package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
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
    protected Parameterization parameterization;

    protected double absoluteTolerance;
    protected double relativeTolerance;
    protected double integrationMinStep;
    protected double integrationMaxStep;

    public IntervalODESystem(Parameterization parameterization, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        this.parameterization = parameterization;
        this.intervals = intervals;
        this.integrationMinStep = this.parameterization.getTotalProcessLength() * 1e-15;
        this.integrationMaxStep = this.parameterization.getTotalProcessLength() / 5;
        this.absoluteTolerance = absoluteTolerance;
        this.relativeTolerance = relativeTolerance;
    }

    /**
     * Integrates over the system forward in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialStates             the initial states at the interval starts-
     * @param intervals                 the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateForwards(List<double[]> initialStates, List<Interval> intervals, boolean alwaysStartAtInitialState, boolean parallelize) {
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState && parallelize) {

            ForkJoinPool forkJoinPool = new ForkJoinPool();

            try {
                forkJoinPool.submit(() ->
                    IntStream.range(0, intervals.size()).parallel().forEach(i -> {
                        Interval interval = intervals.get(i);
                        double[] state = initialStates.get(i).clone();

                        ContinuousOutputModel outputModel = new ContinuousOutputModel();

                        // this interval might contain multiple parameterization intervals
                        // we loop through the containing parameterization intervals and
                        // integrate them piece-by-piece

                        double currentStart = interval.start();
                        for (int parameterizationInterval : interval.parameterizationIntervals()) {
                            this.handleParameterizationIntervalBoundaryIfNecessary(currentStart, state);

                            double parameterizationEnd = this.parameterization.getIntervalEndTimes()[parameterizationInterval];
                            outputModel.append(
                                    this.integrate(state, currentStart, parameterizationEnd, interval)
                            );

                            currentStart = parameterizationEnd;
                        }

                        this.handleParameterizationIntervalBoundaryIfNecessary(currentStart, state);

                        outputModel.append(
                                this.integrate(state, currentStart, interval.end(), interval)
                        );

                        outputModels[interval.interval()] = outputModel;
                    })).join();
            } catch (NumberIsTooSmallException | IllegalStateException exception) {
                // shutdown all threads in case of an exception
                // the exception is automatically passed upwards to the caller
                try {
                    forkJoinPool.shutdownNow();
                    forkJoinPool.awaitTermination(10, TimeUnit.SECONDS);
                } catch (InterruptedException e) {
                    // another thread caused an exception in the meantime
                    // we ignore it and only throw the original one
                }

                throw exception;
            } finally {
                forkJoinPool.shutdown();
            }

        } else {

            double[] state = initialStates.get(0).clone();

            for (int i = 0; i < intervals.size(); i++) {
                Interval interval = intervals.get(i);

                if (alwaysStartAtInitialState) {
                    state = initialStates.get(i).clone();
                }

                ContinuousOutputModel outputModel = new ContinuousOutputModel();

                // this interval might contain multiple parameterization intervals
                // we loop through the containing parameterization intervals and
                // integrate them piece-by-piece

                double currentStart = interval.start();
                for (int parameterizationInterval : interval.parameterizationIntervals()) {
                    this.handleParameterizationIntervalBoundaryIfNecessary(currentStart, state);

                    double parameterizationEnd = this.parameterization.getIntervalEndTimes()[parameterizationInterval];
                    outputModel.append(
                            this.integrate(state, currentStart, parameterizationEnd, interval)
                    );

                    currentStart = parameterizationEnd;
                }

                this.handleParameterizationIntervalBoundaryIfNecessary(currentStart, state);

                outputModel.append(
                        this.integrate(state, currentStart, interval.end(), interval)
                );

                outputModels[interval.interval()] = outputModel;
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

        if (alwaysStartAtInitialState && parallelize) {

            ForkJoinPool forkJoinPool = new ForkJoinPool();

            try {
                forkJoinPool.submit(() ->
                    IntStream.range(0, intervals.size()).parallel().forEach(i -> {
                    Interval interval = intervals.get(i);
                    double[] state = initialStates.get(i).clone();

                    ContinuousOutputModel outputModel = new ContinuousOutputModel();

                    // this interval might contain multiple parameterization intervals
                    // we loop through the containing parameterization intervals and
                    // integrate them piece-by-piece

                    double currentEnd = interval.end();
                    for (int j = interval.parameterizationIntervals().size() - 1; j >= 0; j--) {
                        this.handleParameterizationIntervalBoundaryIfNecessary(currentEnd, state);

                        int parameterizationInterval = interval.parameterizationIntervals().get(j);
                        double parameterizationEnd = this.parameterization.getIntervalEndTimes()[parameterizationInterval];
                        outputModel.append(
                                this.integrate(state, currentEnd, parameterizationEnd, interval)
                        );

                        currentEnd = parameterizationEnd;
                    }

                    this.handleParameterizationIntervalBoundaryIfNecessary(currentEnd, state);

                    outputModel.append(
                            this.integrate(state, currentEnd, interval.start(), interval)
                    );

                    outputModels[intervals.size() - interval.interval() - 1] = outputModel;
                })).join();
            } catch (NumberIsTooSmallException | IllegalStateException exception) {
                // shutdown all threads in case of an exception
                // the exception is automatically passed upwards to the caller
                try {
                    forkJoinPool.shutdownNow();
                    forkJoinPool.awaitTermination(10, TimeUnit.SECONDS);
                } catch (InterruptedException e) {
                    // another thread caused an exception in the meantime
                    // we ignore it and only throw the original one
                }

                throw exception;
            } finally {
                forkJoinPool.shutdown();
            }

        } else {

            double[] state = initialStates.get(initialStates.size() - 1).clone();
            for (int i = intervals.size() - 1; i >= 0; i--) {
                Interval interval = intervals.get(i);

                if (alwaysStartAtInitialState) {
                    state = initialStates.get(i).clone();
                }

                ContinuousOutputModel outputModel = new ContinuousOutputModel();

                // this interval might contain multiple parameterization intervals
                // we loop through the containing parameterization intervals and
                // integrate them piece-by-piece

                double currentEnd = interval.end();
                for (int j = interval.parameterizationIntervals().size() - 1; j >= 0; j--) {
                    this.handleParameterizationIntervalBoundaryIfNecessary(currentEnd, state);

                    int parameterizationInterval = interval.parameterizationIntervals().get(j);
                    double parameterizationEnd = this.parameterization.getIntervalEndTimes()[parameterizationInterval];
                    outputModel.append(
                            this.integrate(state, currentEnd, parameterizationEnd, interval)
                    );

                    currentEnd = parameterizationEnd;
                }

                this.handleParameterizationIntervalBoundaryIfNecessary(currentEnd, state);

                outputModel.append(
                        this.integrate(state, currentEnd, interval.start(), interval)
                );

                outputModels[intervals.size() - interval.interval() - 1] = outputModel;
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
        integrator.integrate(this, start, initialState, end, initialState);
        integrator.clearStepHandlers();

        return intervalResult;
    }

    /**
     * This method is called on parameterization boundaries. Inherit this method for custom handling of these boundaries.
     */
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        // do nothing
    }

    protected void handleParameterizationIntervalBoundary(double boundaryTime, double[] state) {
        this.handleParameterizationIntervalBoundary(boundaryTime, this.getCurrentParameterizationInterval(boundaryTime), this.getCurrentParameterizationInterval(boundaryTime), state);
    }

    protected void handleParameterizationIntervalBoundaryIfNecessary(double boundaryTime, double[] state) {
        if (this.isParameterizationIntervalBoundary(boundaryTime)) {
            this.handleParameterizationIntervalBoundary(boundaryTime, state);
        }
    }

    public int getCurrentParameterizationInterval(double time) {
        return this.parameterization.getIntervalIndex(time);
    }

    boolean isParameterizationIntervalBoundary(double time) {
        for (int i = 0; i < this.parameterization.getTotalIntervalCount() - 1; i++) {
            double endTime = this.parameterization.getIntervalEndTimes()[i];
            if (Utils.equalWithPrecision(endTime, time)) return true;
        }
        return false;
    }

    protected IntervalODESystem copy() {
        return this;
    }
}
