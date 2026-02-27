package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;


import java.util.ArrayList;
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

    int seed;
    double maxConditionNumber;

    public FlowODESystem(
            Parameterization parameterization,
            ExtinctionProbabilities extinctionProbabilities,
            List<Interval> intervals,
            double absoluteTolerance,
            double relativeTolerance,
            int seed,
            double maxConditionNumber
    ) {
        super(parameterization, intervals, absoluteTolerance, relativeTolerance);
        this.extinctionProbabilities = extinctionProbabilities;

        this.birthRates = this.parameterization.getBirthRates();
        this.deathRates = this.parameterization.getDeathRates();
        this.samplingRates = this.parameterization.getSamplingRates();
        this.crossBirthRates = this.parameterization.getCrossBirthRates();
        this.migrationRates = this.parameterization.getMigRates();

        this.seed = seed;
        this.maxConditionNumber = maxConditionNumber;

        this.timeInvariantSystemMatrices = new RealMatrix[this.parameterization.getTotalIntervalCount()];

        for (int i = 0; i < this.parameterization.getTotalIntervalCount(); i++) {
            this.timeInvariantSystemMatrices[i] = this.buildTimeInvariantSystemMatrix(i);
        }
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
                    this.deathRates[interval][i] + this.samplingRates[interval][i] + this.birthRates[interval][i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
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
        ContinuousOutputModel extinctionOutputModel = this.extinctionProbabilities.getOutputModel(t);

        synchronized (extinctionOutputModel) {
            double[] extinctProbabilities = this.extinctionProbabilities.unsafeGetProbability(extinctionOutputModel, t);
            int interval = this.getCurrentParameterizationInterval(t);

            for (int i = 0; i < parameterization.getNTypes(); i++) {
                system.addToEntry(
                        i,
                        i,
                        -2 * this.birthRates[interval][i] * extinctProbabilities[i]
                );

                for (int j = 0; j < parameterization.getNTypes(); j++) {
                    system.addToEntry(
                            i,
                            i,
                            -this.crossBirthRates[interval][i][j] * extinctProbabilities[j]
                    );

                    system.addToEntry(
                            i,
                            j,
                            -this.crossBirthRates[interval][i][j] * extinctProbabilities[i]
                    );
                }
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
        if (Double.isNaN(t)) {
            throw new IllegalStateException("NaN detected during integration.");
        }

        int numTypes = this.parameterization.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);
        RealMatrix yDotMatrix = systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            for (int j = 0; j < this.parameterization.getNTypes(); j++) {
                state[i * this.parameterization.getNTypes() + j] *= (1 - this.parameterization.getRhoValues()[newInterval][i]);
            }
        }
    }

    List<double[]> getInitialStates(String initialMatrixStrategy, List<Interval> intervals) {
        return switch (initialMatrixStrategy) {
            case "random" -> {
                List<double[]> arrays = new ArrayList<>();
                for (Interval interval : intervals) {
                    double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                    RealMatrix matrix = Utils.getRandomMatrix(
                            this.parameterization.getNTypes(), this.seed + interval.interval()
                    );
                    Utils.fillArray(matrix, array);
                    arrays.add(array);
                }

                yield arrays;
            }
            case "identity" -> {
                RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());

                double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                Utils.fillArray(matrix, array);

                List<double[]> arrays = new ArrayList<>();
                for (Interval ignored : intervals) {
                    arrays.add(array);
                }

                yield arrays;
            }
            case "average_inverse" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    double h = interval.end() - interval.start();
                    RealMatrix startA = this.buildSystemMatrix(interval.start() + bdmmprime.util.Utils.globalPrecisionThreshold);
                    RealMatrix threeQuarterA = this.buildSystemMatrix(3.0 * (interval.start() + interval.end()) / 4.0);
                    RealMatrix midA = this.buildSystemMatrix((interval.start() + interval.end()) / 2.0);
                    RealMatrix endA = this.buildSystemMatrix(interval.end() - bdmmprime.util.Utils.globalPrecisionThreshold);

                    RealMatrix startInvX = MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());
                    RealMatrix midInvX = Utils.expm(
                            midA.add(threeQuarterA.scalarMultiply(4)).add(endA).scalarMultiply(h / 2.0 / 6.0)
                    );
                    RealMatrix endInvX = Utils.expm(
                            startA.add(midA.scalarMultiply(4)).add(endA).scalarMultiply(h / 6.0)
                    );

                    RealMatrix averageInvX = startInvX.add(midInvX.scalarMultiply(4)).add(endInvX).scalarMultiply(1.0 / 6.0);

                    double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                    Utils.fillArray(averageInvX, array);

                    arrays.add(array);
                }

                yield arrays;
            }
            case "ed" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    RealMatrix startA = this.buildSystemMatrix(interval.start() + bdmmprime.util.Utils.globalPrecisionThreshold);
                    RealMatrix midA = this.buildSystemMatrix((interval.start() + interval.end()) / 2);
                    RealMatrix endA = this.buildSystemMatrix(interval.end() - bdmmprime.util.Utils.globalPrecisionThreshold);
                    double h = interval.end() - interval.start();

                    RealMatrix simpsonMidA = endA.add(midA.scalarMultiply(4)).add(startA).scalarMultiply(1.0 / 6.0);
                    RealMatrix midX = Utils.expm(simpsonMidA);

                    RealMatrix V = new EigenDecomposition(midX).getV();

                    double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                    Utils.fillArray(V, array);

                    arrays.add(array);
                }

                yield arrays;
            }
            default -> throw new RuntimeException(
                    "Error: initial state strategy not known. Valid strategies are 'random' and 'identity'."
            );
        };
    }

    /**
     * Calculates the flow integral using the given intervals.
     *
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
        List<double[]> initialStates = this.getInitialStates(initialMatrixStrategy, intervals);

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStates,
                intervals,
                resetInitialStateAtIntervalsBoundaries,
                parallelize
        );

        return new Flow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries
        );
    }

    /**
     * Splits up the stored intervals if numerical issues are expected.
     *
     * @return the new list of intervals.
     */
    @Override
    public List<Interval> splitUpIntervals() {
        List<Interval> newIntervals = new ArrayList<>();

        double logMaxConditionNumber = Math.log(this.maxConditionNumber);

        int currentOldIntervalIdx = this.intervals.size() - 1;
        double currentIntervalEnd = this.intervals.get(this.intervals.size() - 1).end();

        while (0 <= currentOldIntervalIdx) {
            Interval currentOldInterval = this.intervals.get(currentOldIntervalIdx);

            double minNewIntervalStart = currentOldInterval.start();
            double minNewIntervalMid = 0.5 * (minNewIntervalStart + currentIntervalEnd);

            RealMatrix currentEndSystemMatrix = this.buildSystemMatrix(currentIntervalEnd - bdmmprime.util.Utils.globalPrecisionThreshold);
            RealMatrix currentMidSystemMatrix = this.buildSystemMatrix(minNewIntervalMid);
            RealMatrix currentStartSystemMatrix = this.buildSystemMatrix(minNewIntervalStart + bdmmprime.util.Utils.globalPrecisionThreshold);

            RealMatrix simpsonSystemApproximation = currentEndSystemMatrix.add(currentMidSystemMatrix.scalarMultiply(4)).add(currentStartSystemMatrix).scalarMultiply(1.0 / 6.0);

            SingularValueDecomposition decomposition = new SingularValueDecomposition(simpsonSystemApproximation);
            double maxSingularValue = Arrays.stream(decomposition.getSingularValues()).max().orElseThrow();
            double maxIntervalSize = logMaxConditionNumber / (2.0 * maxSingularValue);

            double newIntervalStart = Math.max(currentIntervalEnd - maxIntervalSize, currentOldInterval.start());

            // find containing parameterization interval ends
            List<Integer> containingParameterizationIntervalEnds = new ArrayList<>();
            for (int j = 0; j < parameterization.getTotalIntervalCount(); j++) {
                double parameterizationIntervalEndTime = parameterization.getIntervalEndTimes()[j];
                if (
                        bdmmprime.util.Utils.lessThanWithPrecision(newIntervalStart, parameterizationIntervalEndTime) &&
                                bdmmprime.util.Utils.lessThanWithPrecision(parameterizationIntervalEndTime, currentIntervalEnd)
                ) {
                    containingParameterizationIntervalEnds.add(j);
                }
            }

            Interval newInterval = new Interval(
                newIntervals.size(), containingParameterizationIntervalEnds, newIntervalStart, currentIntervalEnd
            );
            newIntervals.add(0, newInterval);

            if (bdmmprime.util.Utils.equalWithPrecision(newIntervalStart, currentOldInterval.start())) {
                // we reached the end of the current old interval
                // let's go to the next one
                currentOldIntervalIdx--;
            }
            currentIntervalEnd = newIntervalStart;
        }

        // update interval indices

        for (int i = 0; i < newIntervals.size(); i++) {
            Interval interval = newIntervals.get(i);
            Interval intervalWithCorrectIdx = new Interval(
                    i, interval.parameterizationIntervals(), interval.start(), interval.end()
            );
            newIntervals.set(i, intervalWithCorrectIdx);
        }
        this.intervals = newIntervals;

        return this.intervals;
    }
}
