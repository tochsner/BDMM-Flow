package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;


/**
 * This class allows to represent ODE systems where the time is split up into intervals.
 * <p>
 * There are two reasons why a time is split into intervals:
 * <p>
 * -    The parameterization of the BDMM model can specify parameterization intervals, where the ODE
 * boundary conditions change at the parameterization interval boundaries.
 * A concrete implementation of this class can inherit the handleParameterizationIntervalBoundary method
 * to specify these changes.
 * <p>
 * -    For numerical stability, it can be useful to split up the time even more fine-grained and to restart
 * integration in each interval.
 */
public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected List<Interval> intervals;
    protected Set<Integer> parameterizationIntervalBoundaries;
    protected Parameterization parameterization;

    protected double absoluteTolerance;
    protected double relativeTolerance;
    protected double integrationMinStep;
    protected double integrationMaxStep;

    public IntervalODESystem(Parameterization parameterization, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        this.parameterization = parameterization;
        this.intervals = intervals;
        this.parameterizationIntervalBoundaries = IntervalUtils.getParameterizationIntervalBoundaries(
                this.intervals
        );

        integrationMinStep = this.parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep = this.parameterization.getTotalProcessLength() / 5;
        this.absoluteTolerance = absoluteTolerance;
        this.relativeTolerance = relativeTolerance;
    }

    /**
     * Integrates over the system forward in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialState              the initial state at time 0.
     * @param intervals                 the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateForwards(double[] initialState, List<Interval> intervals, boolean alwaysStartAtInitialState, boolean parallelize) {
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState && 8 < intervals.size()) {

            ForkJoinPool pool = new ForkJoinPool();
            try {
                pool.submit(() ->
                        IntStream.range(0, intervals.size()).parallel().forEach(i -> {
                            Interval interval = intervals.get(i);
                            double[] state = initialState.clone();

                            if (this.parameterizationIntervalBoundaries.contains(interval.interval())) {
                                this.handleParameterizationIntervalBoundary(
                                        interval.start(),
                                        interval.parameterizationInterval() - 1,
                                        interval.parameterizationInterval(),
                                        state
                                );
                            }

                            outputModels[interval.interval()] = this.integrate(state, interval.start(), interval.end(), interval);
                        })
                ).join();
            } catch (Exception exception) {
                // shutdown all threads in case of an exception
                // the exception is automatically passed upwards to the caller
                pool.shutdown();
                pool.shutdownNow();

                throw exception;
            }

        } else {

            double[] state = initialState.clone();

            for (Interval interval : intervals) {
                if (alwaysStartAtInitialState) {
                    state = initialState.clone();
                }

                if (this.parameterizationIntervalBoundaries.contains(interval.interval())) {
                    this.handleParameterizationIntervalBoundary(
                            interval.start(),
                            interval.parameterizationInterval() - 1,
                            interval.parameterizationInterval(),
                            state
                    );
                }

                outputModels[interval.interval()] = this.integrate(state, interval.start(), interval.end(), interval);
            }

        }

        return outputModels;
    }

    /**
     * Integrates over the system backwards in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialStates             the initial states at the interval ends.
     * @param intervals                 the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateBackwards(List<double[]> initialStates, List<Interval> intervals, boolean alwaysStartAtInitialState, boolean parallelize) {
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState && 8 < intervals.size()) {

            ForkJoinPool pool = new ForkJoinPool();
            try {
                pool.submit(() ->
                        IntStream.range(0, intervals.size()).parallel().forEach(i -> {
                            Interval interval = intervals.get(i);
                            double[] state = initialStates.get(i).clone();

                            if (this.parameterizationIntervalBoundaries.contains(interval.interval() + 1)) {
                                this.handleParameterizationIntervalBoundary(
                                        interval.end(),
                                        interval.parameterizationInterval() + 1,
                                        interval.parameterizationInterval(),
                                        state
                                );
                            }

                            outputModels[intervals.size() - interval.interval() - 1] = this.integrate(
                                    state, interval.end(), interval.start(), interval
                            );
                        })
                ).join();
            } catch (Exception exception) {
                // shutdown all threads in case of an exception
                // the exception is automatically passed upwards to the caller
                pool.shutdown();
                pool.shutdownNow();

                throw exception;
            }

        } else {

            double[] state = initialStates.get(initialStates.size() - 1).clone();
            for (int i = intervals.size() - 1; i >= 0; i--) {
                Interval interval = intervals.get(i);

                if (alwaysStartAtInitialState) {
                    state = initialStates.get(i).clone();
                }

                if (this.parameterizationIntervalBoundaries.contains(interval.interval() + 1)) {
                    this.handleParameterizationIntervalBoundary(
                            interval.end(),
                            interval.parameterizationInterval() + 1,
                            interval.parameterizationInterval(),
                            state
                    );
                }

                outputModels[intervals.size() - interval.interval() - 1] = this.integrate(
                        state, interval.end(), interval.start(), interval
                );
            }

        }

        return outputModels;
    }

    /**
     * Integrate the system along the given interval from start to end using the given initialState.
     */
    protected ContinuousOutputModel integrate(double[] initialState, double start, double end, Interval interval) {
        ContinuousOutputModel intervalResult = new ContinuousOutputModel();

        DormandPrince853Integrator integrator = new DormandPrince853Integrator(
                this.integrationMinStep, this.integrationMaxStep, this.absoluteTolerance, this.relativeTolerance
        );
        integrator.addStepHandler(intervalResult);
        integrator.integrate(this.constrainToInterval(interval), start, initialState, end, initialState);
        integrator.clearStepHandlers();

        return intervalResult;
    }

    /**
     * This method is called on parameterization boundaries. Inherit this method for custom handling of these boundaries.
     */
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        // do nothing
    }

    public int getCurrentParameterizationInterval(double time) {
        return IntervalUtils.getInterval(time, this.intervals).parameterizationInterval();
    }

    /**
     * Returns the system constraint to the given interval. This is useful when we parallelize
     * over the different intervals to achieve thread-safety.
     */
    abstract protected IntervalODESystem constrainToInterval(Interval interval);
}
