package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
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

    HashMap<Pair<Double, Integer>, RealMatrix> flowCache = new HashMap<>();
    HashMap<Pair<Double, Integer>, DecompositionSolver> decompositionCache = new HashMap<>();
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

        Pair<Double, Integer> startKey = new Pair<>(timeStart, interval);
        DecompositionSolver qr = this.decompositionCache.computeIfAbsent(
                startKey,
                k -> new QRDecomposition(
                        this.getFlow(timeStart, interval), 1e-10
                ).getSolver()
        );

        RealVector likelihoodVectorEnd = new ArrayRealVector(endState);
        RealVector likelihoodVectorStart = qr.solve(
                this.operateFlow(timeEnd, interval, likelihoodVectorEnd)
        );

        return likelihoodVectorStart.toArray();
    }

    /**
     * Calculates the flow at a given time.
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
     * Operates the flow at a given time on the given vector.
     * This method supports when the flow integration was restarted using the same initial state
     * at the beginning of every interval. In this case, the flow is calculated by accumulatively
     * multiplying the end flows of the intervals between startingAtInterval and time.
     *
     * @param time               the time for which to query the flow from.
     * @param startingAtInterval where to start the accumulation of the flow if initial state resetting
     *                           was used.
     * @param vector             the vector to multiply the flow with.
     * @return the flow at the given time multiplied with the vector.
     */
    public RealVector operateFlow(double time, int startingAtInterval, RealVector vector) {
        int timeInterval = this.getInterval(time);

        if (!this.wasInitialStateResetAtEachInterval || startingAtInterval == timeInterval)
            return this.getFlow(this.outputModels[timeInterval], time).operate(vector);

        RealVector accumulatedVector = this.getFlow(this.outputModels[timeInterval], time).operate(vector);

        for (int i = timeInterval - 1; i >= startingAtInterval; i--) {
            RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
            accumulatedVector = this.inverseInitialStates.get(i + 1).operate(accumulatedVector);
            accumulatedVector = flowEnd.operate(accumulatedVector);
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
