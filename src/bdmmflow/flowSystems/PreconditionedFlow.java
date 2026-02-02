package bdmmflow.flowSystems;

import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalUtils;
import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.ArrayList;
import java.util.List;

/**
 * This class is a lightweight wrapper of the result of the Flow ODE integration. It allows to easily query the
 * flow at every point in time and to use the flow to efficiently integrate over a time span.
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class PreconditionedFlow extends Flow {
    RealMatrix[] preconditioners;
    RealMatrix[] inversePreconditioners;
    List<Interval> intervals;

    public PreconditionedFlow(ContinuousOutputModel[] outputModels, int n, List<double[]> initialStateArrays, boolean wasInitialStateResetAtEachInterval, RealMatrix[] preconditioners, RealMatrix[] inversePreconditioners, List<Interval> intervals) {
        super(outputModels, n, initialStateArrays, wasInitialStateResetAtEachInterval);
        this.preconditioners = preconditioners;
        this.inversePreconditioners = inversePreconditioners;
        this.intervals = intervals;
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        Interval interval = IntervalUtils.getInterval(time, this.intervals);
        RealMatrix preconditioner = this.preconditioners[interval.interval()];

        synchronized (output) {
            output.setInterpolatedTime(time);
            return preconditioner.multiply(Utils.toMatrix(output.getInterpolatedState(), this.n));
        }
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
    @Override
    public RealMatrix getFlow(double time, int startingAtInterval) {
        int timeInterval = this.getInterval(time);

        if (!this.wasInitialStateResetAtEachInterval || startingAtInterval == timeInterval)
            return this.getFlow(this.outputModels[timeInterval], time);

        if (this.accumulatedFlowCache[startingAtInterval][timeInterval] == null) {
            RealMatrix accumulatedFlow = MatrixUtils.createRealIdentityMatrix(this.n);

            for (int i = startingAtInterval; i < timeInterval; i++) {
                // get flowEnd in physical space: A(finalTime) = M * Y(finalTime)
                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());

                // convert back to preconditioned space: Y(finalTime) = M^(-1) * A(finalTime)
                Interval interval = IntervalUtils.getInterval(this.outputModels[i].getFinalTime(), this.intervals);
                RealMatrix inversePreconditioner = this.inversePreconditioners[interval.interval()];
                RealMatrix yFlowEnd = inversePreconditioner.multiply(flowEnd);

                // accumulate in preconditioned space: Y(T)^(-1) * Y(finalTime)
                accumulatedFlow = this.inverseInitialStates.get(this.inverseInitialStates.size() - i - 2).multiply(yFlowEnd).multiply(accumulatedFlow);
                this.accumulatedFlowCache[startingAtInterval][i + 1] = accumulatedFlow;
            }
        }

        return (
                this.getFlow(this.outputModels[timeInterval], time)
                        .multiply(this.accumulatedFlowCache[startingAtInterval][timeInterval])
        );
    }

}
