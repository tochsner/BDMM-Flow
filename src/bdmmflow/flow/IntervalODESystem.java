package bdmmflow.flow;

import bdmmflow.flow.intervals.Interval;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected Parameterization param;
    protected FirstOrderIntegrator integrator;
    protected int currentInterval;

    public IntervalODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        this.param = parameterization;

        double integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = this.param.getTotalProcessLength() / 20;

        this.integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
        );
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] state) {
        return this.integrateBackwardsOverIntegrals(state, this.param.getTotalProcessLength(), false);
    }

    public ContinuousOutputModel[] integrateOverIntegrals(double[] initialState, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        this.currentInterval = 0;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        double[] state = initialState.clone();

        for (Interval interval : intervals) {
            if (this.currentInterval != interval.parameterizationInterval()) {
                assert this.currentInterval + 1 == interval.parameterizationInterval();

                this.handleIntervalBoundary(interval.start(), this.currentInterval, this.currentInterval + 1, state);
            }

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            if (alwaysStartAtInitialState) {
                state = initialState.clone();
            }

            assert interval.parameterizationInterval() == this.currentInterval;

            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, interval.start(), state, interval.end(), state);
            integrator.clearStepHandlers();

            outputModels[interval.interval()] = intervalResult;
        }

        return outputModels;
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] initialState, double maxIntervalSize, boolean alwaysStartAtInitialState) {
        this.currentInterval = this.param.getTotalIntervalCount() - 1;

        ArrayList<ContinuousOutputModel> outputModels = new ArrayList<>();

        int currentInterval = this.param.getTotalIntervalCount() - 1;
        double[] state = initialState.clone();

        double endTime = this.param.getIntervalEndTimes()[currentInterval];
        while (true) {
            double intervalStartTime = 0 < currentInterval ? this.param.getIntervalEndTimes()[currentInterval - 1] : 0;

            double timeUntilNextInterval = endTime - intervalStartTime;
            boolean useIntervalStart = timeUntilNextInterval < maxIntervalSize || Utils.equalWithPrecision(maxIntervalSize, timeUntilNextInterval);
            double startTime = useIntervalStart ? intervalStartTime : endTime - maxIntervalSize;

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            if (alwaysStartAtInitialState)
                state = initialState.clone();

            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, endTime, state, startTime, state);
            integrator.clearStepHandlers();

            if (useIntervalStart && 0 < currentInterval) {
                this.handleIntervalBoundary(startTime, currentInterval, currentInterval - 1, state);
            }

            outputModels.add(intervalResult);

            if (useIntervalStart)
                currentInterval--;

            endTime = startTime;

            if (startTime <= 0) break;
        }

        return outputModels.toArray(new ContinuousOutputModel[0]);
    }

    protected void handleIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        this.currentInterval = newInterval;
    }

}
