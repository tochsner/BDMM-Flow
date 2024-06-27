package bdmmflow.flow;

import bdmmflow.flow.intervals.Interval;
import bdmmflow.flow.intervals.IntervalUtils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.List;


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

    public ContinuousOutputModel[] integrateOverIntegrals(double[] state) {
        return this.integrateOverIntegrals(
                state,
                IntervalUtils.getIntervals(this.param, this.param.getTotalProcessLength()),
                false
        );
    }

    public ContinuousOutputModel[] integrateOverIntegrals(double[] initialState, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        this.currentParameterizationInterval = 0;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        double[] state = initialState.clone();

        for (Interval interval : intervals) {
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

            if (alwaysStartAtInitialState) {
                state = initialState.clone();
            }

            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, interval.start(), state, interval.end(), state);
            integrator.clearStepHandlers();

            outputModels[interval.interval()] = intervalResult;
        }

        return outputModels;
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] state) {
        return this.integrateBackwardsOverIntegrals(
                state,
                IntervalUtils.getIntervals(this.param, this.param.getTotalProcessLength()),
                false
        );
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] initialState, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        this.currentParameterizationInterval = this.param.getTotalIntervalCount() - 1;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        double[] state = initialState.clone();

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

            outputModels[interval.interval()] = intervalResult;
        }

        return outputModels;
    }

    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        this.currentParameterizationInterval = newInterval;
    }

}
