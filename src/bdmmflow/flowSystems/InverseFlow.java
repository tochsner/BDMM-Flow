package bdmmflow.flowSystems;

import bdmmflow.utils.LRUCache;
import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * This class is a lightweight wrapper of the result of the Inverse Flow ODE integration. It allows to easily query the
 * inverse flow at every point in time and to use the inverse flow to efficiently integrate over a time span.
 * <p>
 * Compared to the normal Flow, the representation used in this class avoids the QR decomposition for the leaf
 * flows when traversing a tree and integrating over the edges.
 * <p>
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class InverseFlow implements IFlow {
    ContinuousOutputModel[] outputModels;

    List<RealMatrix> inverseInitialStates;
    boolean wasInitialStateResetAtEachInterval;
    int n;

    LRUCache<Pair<Double, Integer>, RealMatrix> flowCache = new LRUCache<>(16);
    LRUCache<Pair<Double, Integer>, DecompositionSolver> decompositionCache = new LRUCache<>(16);
    RealMatrix[][] accumulatedFlowCache;

    public InverseFlow(ContinuousOutputModel[] outputModels, int n, List<double[]> initialStateArrays, boolean useIntervals) {
        this.outputModels = outputModels;

        this.n = n;
        this.wasInitialStateResetAtEachInterval = useIntervals;
        this.accumulatedFlowCache = new RealMatrix[outputModels.length][outputModels.length];

        this.inverseInitialStates = new ArrayList<>();
        for (double[] initialStateArray : initialStateArrays) {
            RealMatrix initialStateMatrix = Utils.toMatrix(initialStateArray, this.n);
            RealMatrix inverseInitialStateMatrix = MatrixUtils.inverse(initialStateMatrix);
            this.inverseInitialStates.add(inverseInitialStateMatrix);
        }
    }

    private InverseFlow(ContinuousOutputModel[] outputModels, int n, boolean wasInitialStateResetAtEachInterval, List<RealMatrix> inverseInitialStates) {
        this.outputModels = outputModels;
        this.n = n;
        this.wasInitialStateResetAtEachInterval = wasInitialStateResetAtEachInterval;
        this.inverseInitialStates = inverseInitialStates;
        this.accumulatedFlowCache = new RealMatrix[outputModels.length][outputModels.length];
    }

    @Override
    public double[] integrateUsingFlow(double timeStart, double timeEnd, double[] endState) {
        int interval = this.getInterval(timeStart);

        RealVector likelihoodVectorEnd = new ArrayRealVector(endState);
        RealVector rightHandSide = this.getFlowVectorProduct(timeEnd, interval, likelihoodVectorEnd);

        Pair<Double, Integer> startKey = new Pair<>(timeStart, interval);
        RealVector likelihoodVectorStart = null;
        try {
            DecompositionSolver qr = this.decompositionCache.computeIfAbsent(
                    startKey,
                    k -> new QRDecomposition(
                            this.flowCache.computeIfAbsent(startKey, k_ -> this.getFlow(timeStart, interval)), 1e-10
                    ).getSolver()
            );
            likelihoodVectorStart = qr.solve(rightHandSide);
        } catch (SingularMatrixException e) {
            // we fall back to SVD in case of nearly-singular matrices
            DecompositionSolver svd = new SingularValueDecomposition(
                    this.flowCache.computeIfAbsent(startKey, k_ -> this.getFlow(timeStart, interval))
            ).getSolver();
            likelihoodVectorStart = svd.solve(rightHandSide);
        }

        return likelihoodVectorStart.toArray();
    }

    /**
     * Calculates the flow at a given time.
     * <p>
     * This method supports when the flow integration was restarted using the same initial state
     * at the beginning of every interval. In this case, the flow is calculated by accumulatively
     * multiplying the end flows of the intervals between startingAtInterval and time.
     *
     * @param time               the time for which to query the flow from.
     * @param startingAtInterval where to start the accumulation of the flow if initial state resetting
     *                           was used.
     * @return the flow at the given time.
     */
    public RealMatrix getFlow(double time, int startingAtInterval) {
        int timeInterval = this.getInterval(time);

        if (!this.wasInitialStateResetAtEachInterval || startingAtInterval == timeInterval)
            return this.getFlow(this.outputModels[timeInterval], time);

        if (this.accumulatedFlowCache[startingAtInterval][timeInterval] == null) {
            RealMatrix accumulatedFlow = MatrixUtils.createRealIdentityMatrix(this.n);

            for (int i = startingAtInterval; i < timeInterval; i++) {
                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
                accumulatedFlow = accumulatedFlow.multiply(flowEnd).multiply(this.inverseInitialStates.get(i + 1));
                this.accumulatedFlowCache[startingAtInterval][i + 1] = accumulatedFlow;
            }
        }

        return (
                this.accumulatedFlowCache[startingAtInterval][timeInterval]
                        .multiply(this.getFlow(this.outputModels[timeInterval], time))
        );
    }

    /**
     * Calculates the flow at a given time.
     * <p>
     * This method supports when the flow integration was restarted using the same initial state
     * at the beginning of every interval. In this case, the flow is calculated by accumulatively
     * multiplying the end flows of the intervals between startingAtInterval and time.
     *
     * @param time               the time for which to query the flow from.
     * @param startingAtInterval where to start the accumulation of the flow if initial state resetting
     *                           was used.
     * @return the flow at the given time.
     */
    public RealVector getFlowVectorProduct(double time, int startingAtInterval, RealVector vector) {
        int timeInterval = this.getInterval(time);

        if (!this.wasInitialStateResetAtEachInterval || startingAtInterval == timeInterval)
            return this.getFlow(this.outputModels[timeInterval], time).operate(vector);

        RealVector accumulatedVector = this.getFlow(this.outputModels[timeInterval], time).operate(vector);

        for (int i = timeInterval - 1; i >= startingAtInterval; i--) {
            RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
            accumulatedVector = flowEnd.operate(this.inverseInitialStates.get(i + 1).operate(accumulatedVector));
        }

        return accumulatedVector;
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        synchronized (output) {
            output.setInterpolatedTime(time);
            return Utils.toMatrix(output.getInterpolatedState(), this.n);
        }
    }

    /**
     * Returns the interval corresponding to the given time.
     *
     * @param time the time to get the interval for.
     * @return the interval.
     */
    public int getInterval(double time) {
        if (time < this.outputModels[0].getInitialTime()) {
            return 0;
        }

        for (int i = 0; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (model.getInitialTime() <= time && time <= model.getFinalTime()) {
                return i;
            }
        }

        return this.outputModels.length - 1;
    }

    public InverseFlow copy() {
        ContinuousOutputModel[] clonedOutputModels = new ContinuousOutputModel[this.outputModels.length];

        for (int i = 0; i < this.outputModels.length; i++) {
            clonedOutputModels[i] = new ContinuousOutputModel();
            clonedOutputModels[i].append(this.outputModels[i]);
        }

        return new InverseFlow(
                clonedOutputModels, this.n, this.wasInitialStateResetAtEachInterval, this.inverseInitialStates
        );
    }
}
