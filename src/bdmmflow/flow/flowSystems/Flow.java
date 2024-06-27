package bdmmflow.flow.flowSystems;

import bdmmflow.flow.Utils;
import bdmmflow.flow.extinctionSystem.ExtinctionProbabilities;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 * This class is a lightweight wrapper of the result of the Flow ODE integration. It allows to easily query the
 * flow at every point in time.
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class Flow extends ExtinctionProbabilities {

    RealMatrix inverseInitialState;
    boolean wasInitialStateResetAtEachInterval;

    int n;

    public Flow(ContinuousOutputModel[] flows, int n) {
        this(flows, n, MatrixUtils.createRealIdentityMatrix(n), false);
    }

    public Flow(ContinuousOutputModel[] flows, int n, RealMatrix inverseInitialState, boolean wasInitialStateResetAtEachInterval) {
        super(flows);

        this.n = n;
        this.inverseInitialState = inverseInitialState;
        this.wasInitialStateResetAtEachInterval = wasInitialStateResetAtEachInterval;
    }

    protected RealMatrix getFlow(ContinuousOutputModel output, double time) {
        return Utils.toMatrix(super.getProbability(output, time), this.n);
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
    public RealMatrix getFlow(double time, int startingAtInterval) {
        RealMatrix accumulatedFlow = null;

        for (int i = startingAtInterval; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                if (accumulatedFlow == null) {
                    return this.getFlow(this.outputModels[i], time);
                } else {
                    return this.getFlow(this.outputModels[i], time).multiply(accumulatedFlow);
                }
            }

            if (this.wasInitialStateResetAtEachInterval) {
                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
                accumulatedFlow = this.inverseInitialState.multiply(flowEnd.multiply(accumulatedFlow));
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }
}
