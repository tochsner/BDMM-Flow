package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.ArrayList;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected Parameterization param;
    protected int interval;

    protected FirstOrderIntegrator integrator;

    public IntervalODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        this.param = parameterization;
        double integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = this.param.getTotalProcessLength() / 20;
        this.integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
        );
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] state) {
        return this.integrateBackwardsOverIntegrals(state, this.param.getTotalProcessLength() / 10, false);
    }

    public ContinuousOutputModel[] integrateOverIntegrals(double[] initialState, double maxIntervalSize, boolean alwaysStartAtInitialState) {
        this.interval = 0;

        ArrayList<ContinuousOutputModel> outputModels = new ArrayList<>();

        int currentInterval = 0;
        double[] state = initialState.clone();

        double startTime = 0;
        while (true) {
            double intervalEndTime = this.param.getIntervalEndTimes()[currentInterval];

            double timeUntilNextInterval = intervalEndTime - startTime;
            boolean useIntervalEnd = timeUntilNextInterval < maxIntervalSize || Utils.equalWithPrecision(maxIntervalSize, timeUntilNextInterval);
            double endTime = useIntervalEnd ? intervalEndTime : startTime + maxIntervalSize;

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            if (alwaysStartAtInitialState)
                state = initialState.clone();

            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, startTime, state, endTime, state);

            integrator.clearStepHandlers();

            if (useIntervalEnd && currentInterval != this.param.getTotalIntervalCount() - 1) {
                this.handleIntervalBoundary(startTime, currentInterval, currentInterval + 1, state);
            }

            outputModels.add(intervalResult);

            if (useIntervalEnd)
                currentInterval++;

            startTime = endTime;

            if (endTime >= this.param.getTotalProcessLength()) break;
        }

        return outputModels.toArray(new ContinuousOutputModel[0]);
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] initialState, double maxIntervalSize, boolean alwaysStartAtInitialState) {
        this.interval = this.param.getTotalIntervalCount() - 1;

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
        this.interval = newInterval;
    }

}
