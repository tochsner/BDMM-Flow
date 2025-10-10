package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.utils.LinearTimeInvMatrixSystem;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

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

    List<double[]> getInitialStates(String initialMatrixStrategy, List<Interval> intervals) {
        return switch (initialMatrixStrategy) {
            case "random" -> {

                List<double[]> arrays = new ArrayList<>();
                for (Interval ignored : intervals) {
                    RealMatrix matrix = Utils.getRandomMatrix(this.param.getNTypes());

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

                    RealMatrix initialStateMatrix = computeRegularMinimizer(initialDerivative);

                    double[] initialStateArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(initialStateMatrix, initialStateArray);

                    arrays.add(initialStateArray);
                }

                yield arrays;
            }
            case "probe_heuristic" -> {
                List<double[]> arrays = new ArrayList<>();

                for (Interval interval : intervals) {
                    double probeStartTime = interval.end() - (interval.end() - interval.start()) / 5;
                    double probeEndTime = interval.end();

                    double[] identityMatrixArray = new double[this.param.getNTypes() * this.param.getNTypes()];
                    Utils.fillArray(MatrixUtils.createRealIdentityMatrix(this.param.getNTypes()), identityMatrixArray);

                    // integrate the probe

                    ContinuousOutputModel probeIntegration = new ContinuousOutputModel();

                    ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(
                            (probeEndTime - probeStartTime) / 10
                    );

                    double[] state = identityMatrixArray.clone();

                    integrator.addStepHandler(probeIntegration);
                    integrator.integrate(this, probeEndTime, state, probeStartTime, state);
                    integrator.clearStepHandlers();

                    probeIntegration.setInterpolatedTime(probeStartTime);
                    RealMatrix integrated = Utils.toMatrix(probeIntegration.getInterpolatedState(), param.getNTypes());

                    // integrate a simple exponential

                    RealMatrix initialDerivative = buildSystemMatrix(probeEndTime);
                    ContinuousOutputModel expIntegration = new LinearTimeInvMatrixSystem(initialDerivative).integrateBackwards(
                            identityMatrixArray, probeStartTime, probeEndTime
                    );

                    expIntegration.setInterpolatedTime(probeStartTime);
                    RealMatrix expIntegrated = Utils.toMatrix(expIntegration.getInterpolatedState(), param.getNTypes());

                    // minimize the error

                    RealMatrix approximationError = integrated.subtract(expIntegrated);
                    RealMatrix initialStateMatrix = computeRegularMinimizer(approximationError);

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

    public static RealMatrix computeRegularMinimizer(RealMatrix A) {
        // Step 1: Compute SVD
        SingularValueDecomposition svd = new SingularValueDecomposition(A);
        RealMatrix V = svd.getV();               // V matrix
        double[] sigma = svd.getSingularValues(); // singular values

        int n = sigma.length;

        // Step 3: Build Sigma^-1 diagonal matrix
        RealMatrix SigmaInv = MatrixUtils.createRealDiagonalMatrix(new double[n]);
        for (int i = 0; i < n; i++) {
            SigmaInv.setEntry(i, i, 1.0 / sigma[i]);
        }

        // Step 4: Compute X* = g * V * Sigma^-1 * V^T
        double d = 1.0; // desired determinant
        double g_d = Math.pow(d * Arrays.stream(sigma).reduce(1.0, (a,b) -> a*b), 1.0/n);
        RealMatrix Xstar = V.multiply(SigmaInv).multiply(V.transpose()).scalarMultiply(g_d);

        return Xstar;
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

        return new Flow(
                rawOutputs,
                this.param.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries
        );
    }
}
