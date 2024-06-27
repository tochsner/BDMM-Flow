package bdmmflow.flow.extinctionSystem;

import bdmmflow.flow.intervals.IntervalODESystem;
import bdmmprime.parameterization.Parameterization;

/**
 * This class represents the ODE system for extinction probabilities.
 */
public class ExtinctionProbabilitiesODESystem extends IntervalODESystem {

    private final double[][] birthRates;
    private final double[][] deathRates;
    private final double[][] samplingRates;
    private final double[][][] crossBirthRates;
    private final double[][][] migrationRates;

    public ExtinctionProbabilitiesODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, absoluteTolerance, relativeTolerance);

        this.birthRates = this.param.getBirthRates();
        this.deathRates = this.param.getDeathRates();
        this.samplingRates = this.param.getSamplingRates();
        this.crossBirthRates = this.param.getCrossBirthRates();
        this.migrationRates = this.param.getMigRates();
    }

    @Override
    public int getDimension() {
        return this.param.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        for (int i = 0; i < this.param.getNTypes(); i++) {
            yDot[i] = (
                    this.birthRates[currentParameterizationInterval][i] * y[i]
                            + this.deathRates[currentParameterizationInterval][i] * y[i]
                            + this.samplingRates[currentParameterizationInterval][i] * y[i]
                            - this.birthRates[currentParameterizationInterval][i] * y[i] * y[i]
                            - this.deathRates[currentParameterizationInterval][i]
            );

            for (int j = 0; j < this.param.getNTypes(); j++) {
                yDot[i] += (
                        this.crossBirthRates[currentParameterizationInterval][i][j] * (y[i] - y[i] * y[j])
                                + this.migrationRates[currentParameterizationInterval][i][j] * (y[i] - y[j])
                );
            }
        }
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int type = 0; type < param.getNTypes(); type++) {
            state[type] *= (1.0 - param.getRhoValues()[newInterval][type]);
        }
    }
}
