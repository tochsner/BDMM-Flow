package bdmmflow.flow.flowSystems;

import bdmmflow.flow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.flow.Utils;
import bdmmflow.flow.intervals.IntervalODESystem;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;

import java.util.LinkedHashMap;
import java.util.Map;

public class InverseFlowODESystem extends IntervalODESystem {
    private final ExtinctionProbabilities extinctionProbabilities;

    private RealMatrix[] timeInvariantSystemMatrices;

    private double[][] birthRates;
    private double[][] deathRates;
    private double[][] samplingRates;
    private double[][][] crossBirthRates;
    private double[][][] migrationRates;

    public InverseFlowODESystem(
            Parameterization parameterization,
            ExtinctionProbabilities extinctionProbabilities,
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

        for (int i = 0; i < this.param.getTotalIntervalCount(); i++) {
            this.timeInvariantSystemMatrices[i] = this.buildTimeInvariantSystemMatrix(i);
        }
    }

    @Override
    public int getDimension() {
        return param.getNTypes() * param.getNTypes();
    }

    /**
     * Builds the time-invariant part of the system matrix for a given interval. This can be reused.
     */
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

    /**
     * Builds the time-varying part of the system matrix for a given interval. This has to be computed for every
     * time step.
     */
    protected void addTimeVaryingSystemMatrix(double t, RealMatrix system) {
        double[] extinctProbabilities = this.extinctionProbabilities.getProbability(t);

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    2 * this.birthRates[currentParameterizationInterval][i] * extinctProbabilities[i]
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.crossBirthRates[currentParameterizationInterval][i][j] * extinctProbabilities[j]
                );

                system.addToEntry(
                        i,
                        j,
                        this.crossBirthRates[currentParameterizationInterval][i][j] * extinctProbabilities[i]
                );
            }
        }
    }

    /**
     * Builds the system matrix for a given time point.
     */
    protected RealMatrix buildSystemMatrix(double t) {
        int interval = this.param.getIntervalIndex(t);
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
            double timeEnd,
            double[] endState,
            InverseFlow flow
    ) {
        int interval = flow.getInterval(timeStart);

        RealMatrix flowMatrixEnd = flow.getFlow(timeEnd, interval);
        RealVector likelihoodVectorEnd = new ArrayRealVector(endState);

        RealMatrix flowMatrixStart = flow.getFlow(timeStart, interval);
        DecompositionSolver qr = new QRDecomposition(flowMatrixStart).getSolver();

        RealVector likelihoodVectorStart = qr.solve(flowMatrixEnd.operate(likelihoodVectorEnd));

        return likelihoodVectorStart.toArray();
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int i = 0; i < this.param.getNTypes(); i++) {
            for (int j = 0; j < this.param.getNTypes(); j++) {
                state[i * this.param.getNTypes() + j] *= (1 - this.param.getRhoValues()[newInterval][i]);
            }
        }
    }
}
