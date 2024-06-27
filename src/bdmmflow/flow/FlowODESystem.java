package bdmmflow.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import java.util.LinkedHashMap;
import java.util.Map;

public class FlowODESystem extends IntervalODESystem {
    private final ContinuousIntervalOutputModel extinctionProbabilities;

    static class LRUCache<K, V> extends LinkedHashMap<K, V> {
        private final int cacheSize;

        public LRUCache(int cacheSize) {
            super(20, 0.75F, true);
            this.cacheSize = cacheSize;
        }

        protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
            return size() >= cacheSize;
        }
    }

    static FlowODESystem.LRUCache<Double, DecompositionSolver> cache = new FlowODESystem.LRUCache<>(20);
    private RealMatrix[] timeInvariantSystemMatrices;

    private double[][] birthRates;
    private double[][] deathRates;
    private double[][] samplingRates;
    private double[][][] crossBirthRates;
    private double[][][] migrationRates;

    public FlowODESystem(
            Parameterization parameterization,
            ContinuousIntervalOutputModel extinctionProbabilities,
            double absoluteTolerance,
            double relativeTolerance
    ) {
        super(parameterization, absoluteTolerance, relativeTolerance);
        this.extinctionProbabilities = extinctionProbabilities;

        this.birthRates = this.param.getBirthRates();
        this.deathRates = this.param.getDeathRates();
        this.samplingRates = this.param.getSamplingRates();
        this.crossBirthRates = this.param.getCrossBirthRates();
        this.migrationRates = this.param.getMigRates();

        this.timeInvariantSystemMatrices = new RealMatrix[this.param.getTotalIntervalCount()];

        cache.clear();
    }

    @Override
    public int getDimension() {
        return param.getNTypes() * param.getNTypes();
    }

    protected RealMatrix buildTimeInvariantSystemMatrix(int interval) {
        RealMatrix system = new BlockRealMatrix(param.getNTypes(), param.getNTypes());

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    -this.deathRates[interval][i] - this.samplingRates[interval][i] - this.birthRates[interval][i]
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        -this.migrationRates[interval][i][j] - this.crossBirthRates[interval][i][j]
                );
                system.addToEntry(
                        i,
                        j,
                        this.migrationRates[interval][i][j]
                );
            }
        }

        return system;
    }

    protected void addTimeVaryingSystemMatrix(double t, RealMatrix system) {
        double[] extinctProbabilities = this.extinctionProbabilities.getOutput(t);

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    2 * this.birthRates[currentInterval][i] * extinctProbabilities[i]
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.crossBirthRates[currentInterval][i][j] * extinctProbabilities[j]
                );

                system.addToEntry(
                        i,
                        j,
                        this.crossBirthRates[currentInterval][i][j] * extinctProbabilities[i]
                );
            }
        }
    }

    protected RealMatrix buildSystemMatrix(double t) {
        int interval = this.param.getIntervalIndex(t);

        if (this.timeInvariantSystemMatrices[interval] == null)
            this.timeInvariantSystemMatrices[interval] = this.buildTimeInvariantSystemMatrix(interval);

        RealMatrix systemMatrix = this.timeInvariantSystemMatrices[interval].copy();
        this.addTimeVaryingSystemMatrix(t, systemMatrix);

        return systemMatrix;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int numTypes = this.param.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);

        RealMatrix yDotMatrix = yMatrix.multiply(systemMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }

    public static double[] integrateUsingFlow(
            double timeStart,
            int intervalStart,
            double timeEnd,
            int intervalEnd,
            double[] initialState,
            Flow flow
    ) {
        int interval = flow.getInterval(timeStart);

        RealMatrix flowMatrixEnd = flow.getFlow(timeEnd, interval, true, false);
        RealVector likelihoodVectorEnd = new ArrayRealVector(initialState, false);

        DecompositionSolver qr;

        if (cache.containsKey(timeStart)) {
            qr = cache.get(timeStart);
        } else {
            RealMatrix flowMatrixStart = flow.getFlow(timeStart, interval, false, true);
            qr = new QRDecomposition(flowMatrixStart).getSolver();
            cache.put(timeStart, qr);
        }

        RealVector likelihoodVectorStart = qr.solve(flowMatrixEnd.operate(likelihoodVectorEnd));

        return likelihoodVectorStart.toArray();
    }

    @Override
    protected void handleIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int i = 0; i < this.param.getNTypes(); i++) {
            for (int j = 0; j < this.param.getNTypes(); j++) {
                state[i * this.param.getNTypes() + j] *= (1 - this.param.getRhoValues()[newInterval][i]);
            }
        }
    }
}
