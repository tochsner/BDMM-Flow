package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;


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

    double absoluteTolerance;
    double relativeTolerance;
    double integrationMinStep;
    double integrationMaxStep;

    public IntervalODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        this.param = parameterization;

        integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        integrationMaxStep = this.param.getTotalProcessLength() / 5;
        this.absoluteTolerance = absoluteTolerance;
        this.relativeTolerance = relativeTolerance;

        this.integrator = new DormandPrince853Integrator(
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
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState) {
            Set<Integer> parameterizationIntervals = new HashSet<>();
            int currentParameterizationInterval = 0;
            for (int i = 0; i < intervals.size(); i++) {
                Interval interval = intervals.get(i);
                if (currentParameterizationInterval != interval.parameterizationInterval()) {
                    assert currentParameterizationInterval + 1 == interval.parameterizationInterval();
                    currentParameterizationInterval = interval.parameterizationInterval();
                    parameterizationIntervals.add(i);
                }
            }

            IntStream.range(0, intervals.size()).parallel().forEach(i -> {
                Interval interval = intervals.get(i);
                double[] state = initialState.clone();

                if (parameterizationIntervals.contains(i)) {
                    this.handleParameterizationIntervalBoundary(
                            interval.start(),
                            interval.parameterizationInterval() - 1,
                            interval.parameterizationInterval(),
                            state
                    );
                }

                ContinuousOutputModel intervalResult = new ContinuousOutputModel();

                DormandPrince853Integrator integrator = new DormandPrince853Integrator(
                        integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
                );
                integrator.addStepHandler(intervalResult);
                integrator.integrate(this, interval.start(), state, interval.end(), state);
                integrator.clearStepHandlers();

                outputModels[interval.interval()] = intervalResult;
            });

        } else {

            this.currentParameterizationInterval = 0;

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
                List.of(initialState),
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
     * @param initialStates the initial states at the interval ends.
     * @param intervals the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateBackwards(List<double[]> initialStates, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState) {
            Set<Integer> parameterizationIntervals = new HashSet<>();
            double currentParameterizationInterval = this.param.getTotalIntervalCount() - 1;
            for (int i = intervals.size() - 1; i >= 0; i--) {
                Interval interval = intervals.get(i);
                if (currentParameterizationInterval != interval.parameterizationInterval()) {
                    assert currentParameterizationInterval - 1 == interval.parameterizationInterval();
                    currentParameterizationInterval = interval.parameterizationInterval();
                    parameterizationIntervals.add(i);
                }
            }

            IntStream.range(0, intervals.size()).parallel().forEach(i -> {
                Interval interval = intervals.get(i);
                double[] state = initialStates.get(i).clone();

                if (parameterizationIntervals.contains(i)) {
                    this.handleParameterizationIntervalBoundary(
                            interval.end(),
                            interval.parameterizationInterval() + 1,
                            interval.parameterizationInterval(),
                            state
                    );
                }

                ContinuousOutputModel intervalResult = new ContinuousOutputModel();

                DormandPrince853Integrator integrator = new DormandPrince853Integrator(
                        integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
                );
                integrator.addStepHandler(intervalResult);
                integrator.integrate(this, interval.end(), state, interval.start(), state);
                integrator.clearStepHandlers();

                outputModels[intervals.size() - interval.interval() - 1] = intervalResult;
            });

        } else {

            this.currentParameterizationInterval = this.param.getTotalIntervalCount() - 1;

            double[] state = initialStates.get(initialStates.size() - 1).clone();
            for (int i = intervals.size() - 1; i >= 0; i--) {
                Interval interval = intervals.get(i);

                if (alwaysStartAtInitialState) {
                    state = initialStates.get(i).clone();
                }

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

                integrator.addStepHandler(intervalResult);
                integrator.integrate(this, interval.end(), state, interval.start(), state);
                integrator.clearStepHandlers();

                outputModels[intervals.size() - interval.interval() - 1] = intervalResult;
            }
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
