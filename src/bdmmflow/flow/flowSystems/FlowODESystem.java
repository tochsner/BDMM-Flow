package bdmmflow.flow.flowSystems;

import bdmmflow.flow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.flow.intervals.Interval;
import bdmmflow.flow.intervals.IntervalODESystem;
import bdmmflow.flow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.Arrays;
import java.util.List;

/**
 * This class represents the classical flow ODE.
 */
public class FlowODESystem extends IntervalODESystem implements IFlowODESystem {
    final ExtinctionProbabilities extinctionProbabilities;

    final RealMatrix[] timeInvariantSystemMatrices;

    final double[][] birthRates;
    final double[][] deathRates;
    final double[][] samplingRates;
    final double[][][] crossBirthRates;
    final double[][][] migrationRates;

    public FlowODESystem(
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
    RealMatrix buildTimeInvariantSystemMatrix(int interval) {
        RealMatrix system = new BlockRealMatrix(param.getNTypes(), param.getNTypes());

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    this.deathRates[interval][i] + this.samplingRates[interval][i] + this.birthRates[interval][i]
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.migrationRates[interval][i][j] + this.crossBirthRates[interval][i][j]
                );
                system.addToEntry(
                        i,
                        j,
                        -this.migrationRates[interval][i][j]
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

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    -2 * this.birthRates[currentParameterizationInterval][i] * extinctProbabilities[i]
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        -this.crossBirthRates[currentParameterizationInterval][i][j] * extinctProbabilities[j]
                );

                system.addToEntry(
                        i,
                        j,
                        -this.crossBirthRates[currentParameterizationInterval][i][j] * extinctProbabilities[i]
                );
            }
        }
    }

    /**
     * Builds the system matrix for a given time point.
     */
    RealMatrix buildSystemMatrix(double t) {
        int interval = this.param.getIntervalIndex(t);
        RealMatrix systemMatrix = this.timeInvariantSystemMatrices[interval].copy();
        this.addTimeVaryingSystemMatrix(t, systemMatrix);

        return systemMatrix;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int numTypes = this.param.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);

        RealMatrix yDotMatrix = systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
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

    /**
     * Calculates the flow integral using the given intervals.
     * @param intervals the intervals to use when integrating over the flow.
     * @return the calculated flow.
     */
    @Override
    public IFlow calculateFlowIntegral(
            List<Interval> intervals,
            RealMatrix initialState,
            RealMatrix inverseInitialState,
            boolean resetInitialStateAtIntervalsBoundaries
    ) {
        double[] initialStateArray = new double[this.param.getNTypes() * this.param.getNTypes()];
        Utils.fillArray(initialState, initialStateArray);

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStateArray,
                intervals,
                resetInitialStateAtIntervalsBoundaries
        );

        return new Flow(
                rawOutputs,
                this.param.getNTypes(),
                inverseInitialState,
                resetInitialStateAtIntervalsBoundaries
        );
    }
}
