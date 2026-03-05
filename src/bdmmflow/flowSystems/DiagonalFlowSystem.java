package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

public class DiagonalFlowSystem extends FlowODESystem {

    public DiagonalFlowSystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance, int seed, double maxConditionNumber) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance, seed, maxConditionNumber, false);
    }

    @Override
    public int getDimension() {
        return parameterization.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        ContinuousOutputModel extinctionOutputModel = this.extinctionProbabilities.getOutputModel(t);

        synchronized (extinctionOutputModel) {
            double[] extinctProbabilities = this.extinctionProbabilities.unsafeGetProbability(extinctionOutputModel, t);
            int interval = this.getCurrentParameterizationInterval(t);

            for (int i = 0; i < parameterization.getNTypes(); i++) {
                yDot[i] = -2 * this.birthRates[interval][i] * extinctProbabilities[i];
            }
        }
    }

}
