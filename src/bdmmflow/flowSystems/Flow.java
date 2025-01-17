package bdmmflow.flowSystems;

import bdmmflow.utils.LRUCache;
import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.util.Pair;

/**
 * This class is a lightweight wrapper of the result of the Flow ODE integration. It allows to easily query the
 * flow at every point in time and to use the flow to efficiently integrate over a time span.
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class Flow implements IFlow {
    ContinuousOutputModel[] outputModels;

    RealMatrix inverseInitialState;
    boolean wasInitialStateResetAtEachInterval;
    int n;

    LRUCache<Pair<Double, Integer>, RealMatrix> cache = new LRUCache<>(16);
    RealMatrix[][] accumulatedFlowCache;

    public Flow(ContinuousOutputModel[] flows, int n) {
        this(flows, n, MatrixUtils.createRealIdentityMatrix(n), false);
    }

    public Flow(ContinuousOutputModel[] outputModels, int n, RealMatrix inverseInitialState, boolean wasInitialStateResetAtEachInterval) {
        this.outputModels = outputModels;
        this.n = n;
        this.inverseInitialState = inverseInitialState;
        this.wasInitialStateResetAtEachInterval = wasInitialStateResetAtEachInterval;
        this.accumulatedFlowCache = new RealMatrix[outputModels.length][outputModels.length];
    }

    /**
     * Allows to integrate over an edge of a tree using the pre-computed flow.
     *
     * @param timeStart the time of the node closer to the root.
     * @param timeEnd   the time of the node closer to the leaves.
     * @param endState  the initial state at the node closer to the leaves.
     * @return the integration result at the time of the node closer to the root.
     */
    @Override
    public double[] integrateUsingFlow(double timeStart, double timeEnd, double[] endState) {
        int intervalEnd = this.getInterval(timeEnd);

        Pair<Double, Integer> startKey = new Pair<>(timeStart, intervalEnd);
        Pair<Double, Integer> endKey = new Pair<>(timeEnd, intervalEnd);

        RealMatrix flowMatrixStart = this.cache.computeIfAbsent(startKey, k -> this.getFlow(timeStart, intervalEnd));
        RealMatrix flowMatrixEnd = this.cache.computeIfAbsent(endKey, k -> this.getFlow(timeEnd, intervalEnd));

        RealVector likelihoodVectorEnd = Utils.toVector(endState);

        DecompositionSolver linearSolver = new QRDecomposition(flowMatrixEnd).getSolver();
        RealVector solution = linearSolver.solve(likelihoodVectorEnd);

        RealVector likelihoodVectorStart = flowMatrixStart.operate(solution);

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
                accumulatedFlow = this.inverseInitialState.multiply(flowEnd.multiply(accumulatedFlow));
            }

            this.accumulatedFlowCache[startingAtInterval][timeInterval] = accumulatedFlow;
        }

        return (
                this.getFlow(this.outputModels[timeInterval], time)
                        .multiply(this.accumulatedFlowCache[startingAtInterval][timeInterval])
        );
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return Utils.toMatrix(output.getInterpolatedState(), this.n);
    }

    /**
     * Returns the interval corresponding to the given time.
     *
     * @param time the time to get the interval for.
     * @return the interval.
     */
    public int getInterval(double time) {
        if (this.outputModels[0].getInitialTime() < time) {
            return 0;
        }

        for (int i = 0; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (model.getFinalTime() < time && time <= model.getInitialTime()) {
                return i;
            }
        }

        return this.outputModels.length - 1;
    }
}
