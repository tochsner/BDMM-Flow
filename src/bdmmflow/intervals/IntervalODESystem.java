package bdmmflow.intervals;

import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalUtils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.List;


/**
 * This class allows to represent ODE systems where the time is split up into intervals.
 * <p>
 * There are two reasons why a time is split into intervals:
 * <p>
 * -    The parameterization of the BDMM model can specify parameterization intervals, where the ODE
 *      boundary conditions change at the parameterization interval boundaries.
 *      A concrete implementation of this class can inherit the handleParameterizationIntervalBoundary method
 *      to specify these changes.
 * <p>
 * -    For numerical stability, it can be useful to split up the time even more fine-grained and to restart
 *      integration in each interval.
 */
public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected int currentParameterizationInterval;
    protected Parameterization param;
    protected FirstOrderIntegrator integrator;

    public IntervalODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        this.param = parameterization;

        double integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = this.param.getTotalProcessLength() / 20;

        this.integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
        );
    }

    /**
     * Integrates over the system forward in time. The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialState the initial state at time 0.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateForwards(double[] initialState) {
        return this.integrateForwards(
                initialState,
                IntervalUtils.getIntervals(this.param, this.param.getTotalProcessLength()),
                false
        );
    }

    /**
     * Integrates over the system forward in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialState the initial state at time 0.
     * @param intervals the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateForwards(double[] initialState, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        this.currentParameterizationInterval = 0;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        double[] state = initialState.clone();

        for (Interval interval : intervals) {
            if (alwaysStartAtInitialState) {
                state = initialState.clone();
            }

            if (this.currentParameterizationInterval != interval.parameterizationInterval()) {
                assert this.currentParameterizationInterval + 1 == interval.parameterizationInterval();

                this.handleParameterizationIntervalBoundary(
                        interval.start(),
                        this.currentParameterizationInterval,
                        this.currentParameterizationInterval + 1,
                        state
                );
            }

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, interval.start(), state, interval.end(), state);
            integrator.clearStepHandlers();

            outputModels[interval.interval()] = intervalResult;
        }

        return outputModels;
    }

    /**
     * Integrates over the system backwards in time. The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialState the initial state at time totalProcessLength.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateBackwards(double[] initialState) {
        return this.integrateBackwards(
                initialState,
                IntervalUtils.getIntervals(this.param, this.param.getTotalProcessLength()),
                false
        );
    }

    /**
     * Integrates over the system backwards in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialState the initial state at time totalProcessLength.
     * @param intervals the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateBackwards(double[] initialState, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        this.currentParameterizationInterval = this.param.getTotalIntervalCount() - 1;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        double[] state = initialState.clone();

        this.handleParameterizationIntervalBoundary(
                this.param.getTotalProcessLength(),
                this.currentParameterizationInterval + 1,
                this.currentParameterizationInterval,
                state
        );

        for (int i = intervals.size() - 1; i >= 0; i--) {
            Interval interval = intervals.get(i);

            if (this.currentParameterizationInterval != interval.parameterizationInterval()) {
                assert this.currentParameterizationInterval - 1 == interval.parameterizationInterval();

                this.handleParameterizationIntervalBoundary(
                        interval.end(),
                        this.currentParameterizationInterval,
                        this.currentParameterizationInterval - 1,
                        state
                );
            }

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            if (alwaysStartAtInitialState) {
                state = initialState.clone();
            }

            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, interval.end(), state, interval.start(), state);
            integrator.clearStepHandlers();

            outputModels[intervals.size() - interval.interval() - 1] = intervalResult;
        }

        return outputModels;
    }

    /**
     * This method is called on parameterization boundaries. Inherit this method for custom handling of these boundaries.
     */
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        this.currentParameterizationInterval = newInterval;
    }

}
