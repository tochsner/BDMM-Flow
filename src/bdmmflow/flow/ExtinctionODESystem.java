package bdmmflow.flow;

import bdmmprime.parameterization.Parameterization;


public class ExtinctionODESystem extends IntervalODESystem {

    private final double[][] birthRates;
    private final double[][] deathRates;
    private final double[][] samplingRates;
    private final double[][][] crossBirthRates;
    private final double[][][] migrationRates;

    public ExtinctionODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
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
            yDot[i] = this.birthRates[currentInterval][i] * y[i];
            yDot[i] += this.deathRates[currentInterval][i] * y[i];
            yDot[i] += this.samplingRates[currentInterval][i] * y[i];

            yDot[i] -= this.birthRates[currentInterval][i] * y[i] * y[i];
            yDot[i] -= this.deathRates[currentInterval][i];

            for (int j = 0; j < this.param.getNTypes(); j++) {
                yDot[i] += this.crossBirthRates[currentInterval][i][j] * (y[i] - y[i] * y[j]);
                yDot[i] += this.migrationRates[currentInterval][i][j] * (y[i] - y[j]);
            }
        }
    }

    @Override
    protected void handleIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int type = 0; type < param.getNTypes(); type++) {
            state[type] *= (1.0 - param.getRhoValues()[newInterval][type]);
        }
    }
}
