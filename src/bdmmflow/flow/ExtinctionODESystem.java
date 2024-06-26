package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;


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
            yDot[i] = this.birthRates[interval][i] * y[i];
            yDot[i] += this.deathRates[interval][i] * y[i];
            yDot[i] += this.samplingRates[interval][i] * y[i];

            yDot[i] -= this.birthRates[interval][i] * y[i] * y[i];
            yDot[i] -= this.deathRates[interval][i];

            for (int j = 0; j < this.param.getNTypes(); j++) {
                yDot[i] += this.crossBirthRates[interval][i][j] * (y[i] - y[i] * y[j]);
                yDot[i] += this.migrationRates[interval][i][j] * (y[i] - y[j]);
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
