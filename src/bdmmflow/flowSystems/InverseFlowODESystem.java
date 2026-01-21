package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.utils.LRUCache;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

import static bdmmflow.utils.Utils.*;

/**
 * This class represents the ODE that has the inverse flow as a solution.
 */
public class InverseFlowODESystem extends IntervalODESystem implements IFlowODESystem {
    final ExtinctionProbabilities extinctionProbabilities;

    final RealMatrix[] timeInvariantSystemMatrices;

    final double[][] birthRates;
    final double[][] deathRates;
    final double[][] samplingRates;
    final double[][][] crossBirthRates;
    final double[][][] migrationRates;

    int seed;

    static LRUCache<Double, DecompositionSolver> cache = new LRUCache<>(16);

    public InverseFlowODESystem(
            Parameterization parameterization,
            ExtinctionProbabilities extinctionProbabilities,
            List<Interval> intervals,
            double absoluteTolerance,
            double relativeTolerance,
            int seed
    ) {
        super(parameterization, intervals, absoluteTolerance, relativeTolerance);
        this.extinctionProbabilities = extinctionProbabilities;

        this.birthRates = this.parameterization.getBirthRates();
        this.deathRates = this.parameterization.getDeathRates();
        this.samplingRates = this.parameterization.getSamplingRates();
        this.crossBirthRates = this.parameterization.getCrossBirthRates();
        this.migrationRates = this.parameterization.getMigRates();

        this.seed = seed;

        this.timeInvariantSystemMatrices = new RealMatrix[this.parameterization.getTotalIntervalCount()];

        for (int i = 0; i < this.parameterization.getTotalIntervalCount(); i++) {
            this.timeInvariantSystemMatrices[i] = this.buildTimeInvariantSystemMatrix(i);
        }

        cache.clear();
    }

    @Override
    public int getDimension() {
        return parameterization.getNTypes() * parameterization.getNTypes();
    }

    /**
     * Builds the time-invariant part of the system matrix for a given interval. This can be reused.
     */
    RealMatrix buildTimeInvariantSystemMatrix(int interval) {
        RealMatrix system = new BlockRealMatrix(parameterization.getNTypes(), parameterization.getNTypes());

        for (int i = 0; i < parameterization.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    -this.deathRates[interval][i] - this.samplingRates[interval][i] - this.birthRates[interval][i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
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
    void addTimeVaryingSystemMatrix(double t, RealMatrix system) {
        double[] extinctProbabilities = this.extinctionProbabilities.getProbability(t);
        int interval = getCurrentParameterizationInterval(t);

        for (int i = 0; i < parameterization.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    2 * this.birthRates[interval][i] * extinctProbabilities[i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.crossBirthRates[interval][i][j] * extinctProbabilities[j]
                );

                system.addToEntry(
                        i,
                        j,
                        this.crossBirthRates[interval][i][j] * extinctProbabilities[i]
                );
            }
        }
    }

    /**
     * Builds the system matrix for a given time point.
     */
    RealMatrix buildSystemMatrix(double t) {
        int interval = this.parameterization.getIntervalIndex(t);
        RealMatrix systemMatrix = this.timeInvariantSystemMatrices[interval].copy();
        this.addTimeVaryingSystemMatrix(t, systemMatrix);

        return systemMatrix;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int numTypes = this.parameterization.getNTypes();

        RealMatrix yMatrix = toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);

        RealMatrix yDotMatrix = yMatrix.multiply(systemMatrix);
        fillArray(yDotMatrix, yDot);
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            for (int j = 0; j < this.parameterization.getNTypes(); j++) {
                state[i * this.parameterization.getNTypes() + j] *= (1 - this.parameterization.getRhoValues()[oldInterval][i]);
            }
        }
    }

    RealMatrix getInitialState(String initialMatrixStrategy) {
        return switch (initialMatrixStrategy) {
            case "random" -> Utils.getRandomMatrix(this.parameterization.getNTypes(), this.seed);
            case "heuristic" -> MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());
            default -> MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());
        };
    }

    /**
     * Calculates the flow integral using the given intervals.
     * @param intervals the intervals to use when integrating over the flow.
     * @return the calculated flow.
     */
    @Override
    public IFlow calculateFlowIntegral(
            List<Interval> intervals,
            String initialMatrixStrategy,
            boolean resetInitialStateAtIntervalsBoundaries,
            boolean parallelize
    ) {
        RealMatrix initialState = this.getInitialState(initialMatrixStrategy);
        double[] initialStateArray = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
        Utils.fillArray(initialState, initialStateArray);

        RealMatrix  inverseInitialState = MatrixUtils.inverse(initialState);

        ContinuousOutputModel[] rawOutputs = this.integrateForwards(
                initialStateArray,
                intervals,
                resetInitialStateAtIntervalsBoundaries,
                parallelize
        );
        return new InverseFlow(
                rawOutputs,
                this.parameterization.getNTypes(),
                inverseInitialState,
                resetInitialStateAtIntervalsBoundaries
        );
    }

    @Override
    public List<Interval> splitUpIntervals() {
        return this.intervals;
    }
}
