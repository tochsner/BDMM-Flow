package bdmmflow.flow.flowSystems;

import bdmmflow.flow.utils.LRUCache;
import bdmmflow.flow.utils.Utils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.util.Pair;

/**
 * This class is a lightweight wrapper of the result of the Flow ODE integration. It allows to easily query the
 * flow at every point in time.
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class InverseFlow implements IFlow {
    ContinuousOutputModel[] outputModels;

    RealMatrix inverseInitialState;
    boolean wasInitialStateResetAtEachInterval;
    int n;

    LRUCache<Pair<Double, Integer>, RealMatrix> flowCache = new LRUCache<>(16);
    LRUCache<Pair<Double, Integer>, DecompositionSolver> decompositionCache = new LRUCache<>(16);

    public InverseFlow(ContinuousOutputModel[] outputModels, int n, RealMatrix inverseInitialState, boolean useIntervals) {
        this.outputModels = outputModels;

        this.n = n;
        this.inverseInitialState = inverseInitialState;
        this.wasInitialStateResetAtEachInterval = useIntervals;
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return Utils.toMatrix(output.getInterpolatedState(), this.n);
    }

    @Override
    public double[] integrateUsingFlow(double timeStart, double timeEnd, double[] endState) {
        int interval = this.getInterval(timeStart);

        Pair<Double, Integer> endKey = new Pair<>(timeEnd, interval);
        RealMatrix flowMatrixEnd = this.flowCache.computeIfAbsent(
                endKey,
                k -> this.getFlow(timeEnd, interval)
        );
        RealVector likelihoodVectorEnd = new ArrayRealVector(endState);

        Pair<Double, Integer> startKey = new Pair<>(timeStart, interval);
        DecompositionSolver qr = this.decompositionCache.computeIfAbsent(
            startKey,
            k -> new QRDecomposition(
                    this.flowCache.computeIfAbsent(startKey, k_ -> this.getFlow(timeStart, interval))
            ).getSolver()
        );

        RealVector likelihoodVectorStart = qr.solve(flowMatrixEnd.operate(likelihoodVectorEnd));

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
        RealMatrix accumulatedFlow = null;

        for (int i = startingAtInterval; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                RealMatrix result;

                if (accumulatedFlow == null) {
                    result = this.getFlow(this.outputModels[i], time);
                } else {
                    result = accumulatedFlow.multiply(this.getFlow(this.outputModels[i], time));
                }

                return result;
            }

            if (this.wasInitialStateResetAtEachInterval) {
                if (accumulatedFlow == null) {
                    accumulatedFlow = MatrixUtils.createRealIdentityMatrix(this.n);
                }

                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
                accumulatedFlow = accumulatedFlow.multiply(flowEnd).multiply(this.inverseInitialState);
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }

    /**
     * Returns the interval corresponding to the given time.
     *
     * @param time the time to get the interval for.
     * @return the interval.
     */
    public int getInterval(double time) {
        for (int i = 0; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                return i;
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }
}
