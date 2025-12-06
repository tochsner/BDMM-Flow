package bdmmflow.flowSystems;

import bdmmflow.benchmark.BenchmarkRun;
import bdmmflow.benchmark.LoggedMetric;
import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.utils.LinearTimeInvMatrixSystem;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.jblas.MatrixFunctions;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.lang.reflect.Field;
import java.util.stream.IntStream;

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

    // this.deathRates[interval][i] + this.samplingRates[interval][i] + this.birthRates[interval][i] ( 1 - 2 * extinctProbabilities[i]) + sum_j( this.crossBirthRates[interval][i][j] (1.0 - 2* extinctProbabilities[j] )


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

    List<double[]> getInitialStates(String initialMatrixStrategy, List<Interval> intervals) {
        return switch (initialMatrixStrategy) {
            case "random" -> {
                RealMatrix matrix = Utils.getRandomMatrix(this.param.getNTypes());

                List<double[]> arrays = new ArrayList<>();
                for (Interval ignored : intervals) {
                    double[] array = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(matrix, array);
                    arrays.add(array);
                }

                yield arrays;
            }
            case "error_heuristic" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    RealMatrix initialDerivative = buildSystemMatrix(interval.end());
                    RealMatrix initialStateMatrix = Utils.computeRegularMinimizer(initialDerivative);

                    double[] initialStateArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(initialStateMatrix, initialStateArray);

                    arrays.add(initialStateArray);
                }

                yield arrays;
            }
            case "taylor_heuristic" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    double probeStartTime = interval.end() - (interval.end() - interval.start()) / 3;
                    double probeEndTime = interval.end();

                    double[] identityMatrixArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    RealMatrix identityMatrix = MatrixUtils.createRealIdentityMatrix(this.param.getNTypes());
                    Utils.fillArray(identityMatrix, identityMatrixArray);

                    // integrate the probe

                    ContinuousOutputModel probeIntegration = new ContinuousOutputModel();

                    EulerIntegrator integrator = new EulerIntegrator(
                            (probeEndTime - probeStartTime) / 5
                    );

                    double[] state = identityMatrixArray.clone();

                    integrator.addStepHandler(probeIntegration);
                    integrator.integrate(this, probeEndTime, state, probeStartTime, state);
                    integrator.clearStepHandlers();

                    probeIntegration.setInterpolatedTime(probeStartTime);
                    RealMatrix integrated = Utils.toMatrix(probeIntegration.getInterpolatedState(), param.getNTypes());

                    // perform taylor approximation

                    RealMatrix initialDerivative = buildSystemMatrix(interval.end());
                    RealMatrix taylorApproximation = identityMatrix.add(initialDerivative.scalarMultiply(probeStartTime - probeEndTime));

                    // minimize the error

                    RealMatrix approximationError = integrated.subtract(taylorApproximation);
                    RealMatrix initialStateMatrix = Utils.computeRegularMinimizer(approximationError);

                    double[] initialStateArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(initialStateMatrix, initialStateArray);

                    arrays.add(initialStateArray);
                }

                yield arrays;
            }
            case "taylor_exp_heuristic" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    double probeStartTime = interval.end() - (interval.end() - interval.start()) / 3;
                    double probeEndTime = interval.end();

                    double[] identityMatrixArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    RealMatrix identityMatrix = MatrixUtils.createRealIdentityMatrix(this.param.getNTypes());
                    Utils.fillArray(identityMatrix, identityMatrixArray);

                    // integrate a simple exponential

                    RealMatrix initialDerivative = buildSystemMatrix(probeEndTime);
                    RealMatrix expIntegrated = Utils.toMatrix(
                            MatrixFunctions.expm(Utils.toMatrix(initialDerivative.scalarMultiply( probeStartTime - probeEndTime)))
                    );

                    // perform taylor approximation

                    initialDerivative = buildSystemMatrix(interval.end());
                    RealMatrix taylorApproximation = identityMatrix.add(initialDerivative.scalarMultiply(probeStartTime - probeEndTime));

                    // minimize the error

                    RealMatrix approximationError = expIntegrated.subtract(taylorApproximation);
                    RealMatrix initialStateMatrix = Utils.computeRegularMinimizer(approximationError);

                    double[] initialStateArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(initialStateMatrix, initialStateArray);

                    arrays.add(initialStateArray);
                }

                yield arrays;
            }
            case "probe_heuristic" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    double probeStartTime = interval.end() - (interval.end() - interval.start()) / 3;
                    double probeEndTime = interval.end();

                    double[] identityMatrixArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(MatrixUtils.createRealIdentityMatrix(this.param.getNTypes()), identityMatrixArray);

                    // integrate the probe

                    ContinuousOutputModel probeIntegration = new ContinuousOutputModel();

                    EulerIntegrator integrator = new EulerIntegrator(
                            (probeEndTime - probeStartTime) / 5
                    );

                    double[] state = identityMatrixArray.clone();

                    integrator.addStepHandler(probeIntegration);
                    integrator.integrate(this, probeEndTime, state, probeStartTime, state);
                    integrator.clearStepHandlers();

                    probeIntegration.setInterpolatedTime(probeStartTime);
                    RealMatrix integrated = Utils.toMatrix(probeIntegration.getInterpolatedState(), param.getNTypes());

                    // integrate a simple exponential

                    RealMatrix initialDerivative = buildSystemMatrix(probeEndTime);
                    RealMatrix expIntegrated = Utils.toMatrix(
                            MatrixFunctions.expm(Utils.toMatrix(initialDerivative.scalarMultiply( probeStartTime - probeEndTime)))
                    );

                    // minimize the error

                    RealMatrix approximationError = expIntegrated.subtract(integrated);
                    RealMatrix initialStateMatrix = Utils.computeRegularMinimizer(approximationError);

                    double[] initialStateArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(initialStateMatrix, initialStateArray);

                    arrays.add(initialStateArray);
                }

                yield arrays;
            }
            default -> {
                RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(this.param.getNTypes());

                double[] array = new double[this.param.getNTypes() * this.param.getNTypes()];
                Utils.fillArray(matrix, array);

                List<double[]> arrays = new ArrayList<>();
                for (Interval ignored : intervals) {
                    arrays.add(array);
                }

                yield arrays;
            }
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
            boolean resetInitialStateAtIntervalsBoundaries
    ) {
        List<double[]> initialStates = this.getInitialStates(initialMatrixStrategy, intervals);

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStates,
                intervals,
                resetInitialStateAtIntervalsBoundaries
        );

        // log num integration steps taken

//        int numIntegrationSteps = 0;
//
//        for (ContinuousOutputModel integrated : rawOutputs) {
//            try {
//                Field field = integrated.getClass().getDeclaredField("index");
//                field.setAccessible(true);
//                numIntegrationSteps += (Integer) field.get(integrated);
//            } catch (Exception e) {
//                throw new RuntimeException(e);
//            }
//        }
//
//        BenchmarkRun.addToMetric("num_steps", Double.valueOf(numIntegrationSteps));

        Flow flow = new Flow(
                rawOutputs,
                this.param.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries
        );

        return flow;
    }
}
