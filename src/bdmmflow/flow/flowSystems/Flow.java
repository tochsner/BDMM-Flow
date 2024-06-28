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
public class Flow implements IFlow {
    ContinuousOutputModel[] outputModels;

    RealMatrix inverseInitialState;
    boolean wasInitialStateResetAtEachInterval;
    int n;

    LRUCache<Pair<Double, Integer>, RealMatrix> cache = new LRUCache<>(16);

    public Flow(ContinuousOutputModel[] flows, int n) {
        this(flows, n, MatrixUtils.createRealIdentityMatrix(n), false);
    }

    public Flow(ContinuousOutputModel[] outputModels, int n, RealMatrix inverseInitialState, boolean wasInitialStateResetAtEachInterval) {
        this.outputModels = outputModels;
        this.n = n;
        this.inverseInitialState = inverseInitialState;
        this.wasInitialStateResetAtEachInterval = wasInitialStateResetAtEachInterval;
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return Utils.toMatrix(output.getInterpolatedState(), this.n);
    }

    /**
     * Returns the interval corresponding to the given time.
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

    /**
     * Calculates the flow at a given time.
     * <p>
     * This method supports when the flow integration was restarted using the same initial state
     * at the beginning of every interval. In this case, the flow is calculated by accumulatively
     * multiplying the end flows of the intervals between startingAtInterval and time.
     * @param time the time for which to query the flow from.
     * @param startingAtInterval where to start the accumulation of the flow if initial state resetting
     *                           was used.
     * @return the flow at the given time.
     */
    @Override
    public RealMatrix getFlow(double time, int startingAtInterval) {
        return this.getFlow(time, startingAtInterval, false, false);
    }

    /**
     * Calculates the flow at a given time.
     * <p>
     * This method supports when the flow integration was restarted using the same initial state
     * at the beginning of every interval. In this case, the flow is calculated by accumulatively
     * multiplying the end flows of the intervals between startingAtInterval and time.
     * @param time the time for which to query the flow from.
     * @param startingAtInterval where to start the accumulation of the flow if initial state resetting
     *                           was used.
     * @return the flow at the given time.
     */
    public RealMatrix getFlow(double time, int startingAtInterval, boolean checkCache, boolean storeInCache) {
        Pair<Double, Integer> key = new Pair<>(time, startingAtInterval);

        if (checkCache && this.cache.containsKey(key))
            return this.cache.get(key);

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
                    result = this.getFlow(this.outputModels[i], time).multiply(accumulatedFlow);
                }

                if (storeInCache) this.cache.put(key, result);

                return result;
            }

            if (this.wasInitialStateResetAtEachInterval) {
                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
                accumulatedFlow = this.inverseInitialState.multiply(flowEnd.multiply(accumulatedFlow));
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }

    /**
     * Allows to integrate over an edge of a tree using the pre-computed flow.
     *
     * @param timeStart the time of the node closer to the root.
     * @param timeEnd  the time of the node closer to the leaves.
     * @param endState the initial state at the node closer to the leaves.
     * @return the integration result at the time of the node closer to the root.
     */
    @Override
    public double[] integrateUsingFlow(double timeStart, double timeEnd, double[] endState) {
        int intervalStart = this.getInterval(timeStart);

        RealMatrix flowMatrixStart = this.getFlow(timeStart, intervalStart, true, true);
        RealMatrix flowMatrixEnd = this.getFlow(timeEnd, intervalStart, true, false);

        RealVector likelihoodVectorEnd = bdmmprime.flow.Utils.toVector(endState);

        DecompositionSolver linearSolver = new QRDecomposition(flowMatrixEnd).getSolver();
        RealVector solution = linearSolver.solve(likelihoodVectorEnd);

        RealVector likelihoodVectorStart = flowMatrixStart.operate(solution);

        return likelihoodVectorStart.toArray();
    }
}
