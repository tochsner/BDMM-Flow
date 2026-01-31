package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.ArrayList;
import java.util.List;

/**
 * This class represents the classical flow ODE.
 */
public class IPFlowODESystem extends FlowODESystem {

    public IPFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance, int seed, double maxConditionNumber) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance, seed, maxConditionNumber);
    }

    /**
     * Builds the system matrix for a given time point.
     */
    RealMatrix buildSystemMatrix(double t, double[] y) {
        int n = this.parameterization.getNTypes();
        int interval = this.parameterization.getIntervalIndex(t);

        RealMatrix system = this.timeInvariantSystemMatrices[interval].copy();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double factor = Math.exp(y[n * n + j] - y[n * n + i]);
                system.multiplyEntry(i, j, factor);
            }
        }

        return system;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        // fill matrix

        int numTypes = this.parameterization.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t, y);
        RealMatrix yDotMatrix = systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);

        // fill diagonal parts

        int n = this.parameterization.getNTypes();

        ContinuousOutputModel extinctionOutputModel = this.extinctionProbabilities.getOutputModel(t);

        synchronized (extinctionOutputModel) {
            double[] extinctProbabilities = this.extinctionProbabilities.unsafeGetProbability(extinctionOutputModel, t);
            int interval = this.getCurrentParameterizationInterval(t);

            for (int i = 0; i < n; i++) {
                yDot[n * n + i] = -2 * this.birthRates[interval][i] * extinctProbabilities[i];
            }
        }
    }

    List<double[]> getInitialStates(String initialMatrixStrategy, List<Interval> intervals) {
        int n = this.parameterization.getNTypes();

        RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(n);

        double[] array = new double[n*n + n];
        Utils.fillArray(matrix, array);

        for (int i = 0; i < n; i++) {
            array[n*n + i] = 1.0;
        }

        List<double[]> arrays = new ArrayList<>();
        for (Interval ignored : intervals) {
            arrays.add(array);
        }

        return arrays;
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

        return new IPFlow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries,
                this.parameterization
        );
    }

    @Override
    public int getDimension() {
        return parameterization.getNTypes() * parameterization.getNTypes() + parameterization.getNTypes();
    }

}
