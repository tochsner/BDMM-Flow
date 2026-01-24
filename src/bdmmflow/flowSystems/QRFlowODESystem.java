package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.ArrayList;
import java.util.List;

/**
 * This class represents the classical flow ODE.
 */
public class QRFlowODESystem extends FlowODESystem {

    public QRFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance, int seed, double maxConditionNumber) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance, seed, maxConditionNumber);
    }

    @Override
    public int getDimension() {
        return 2 * parameterization.getNTypes() * parameterization.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int n = this.parameterization.getNTypes();

        RealMatrix Q = Utils.toMatrix(y, n);
        RealMatrix R = Utils.toMatrix(y, n, n*n);

        RealMatrix A = this.buildSystemMatrix(t);

        // B = Q^T*A*Q
        RealMatrix B = Q.transpose().multiply(A).multiply(Q);

        // add Q_dot = QH to derivative

        RealMatrix H = new BlockRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i > j) {
                    H.setEntry(i, j, B.getEntry(i, j));
                } else if (j > i) {
                    H.setEntry(i, j, -B.getEntry(j, i));
                }
            }
        }

        RealMatrix QDot = Q.multiply(H);
        Utils.fillArray(QDot, yDot);

        // add R_dot = MR to derivative

        RealMatrix M = new BlockRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j > i) {
                    M.setEntry(i, j, B.getEntry(i, j) + B.getEntry(j, i));
                } else if (i == j) {
                    M.setEntry(i, j, B.getEntry(i, j));
                }
            }
        }

        RealMatrix RDot = M.multiply(R);
        Utils.fillArray(RDot, yDot, n*n);
    }

    List<double[]> getInitialStates(List<Interval> intervals) {
        int n = this.parameterization.getNTypes();

        RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());

        double[] array = new double[this.getDimension()];
        Utils.fillArray(matrix, array);
        Utils.fillArray(matrix, array, n*n);

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
        List<double[]> initialStates = this.getInitialStates(intervals);

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStates,
                intervals,
                false,
                parallelize
        );

        return new QRFlow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries
        );
    }

    @Override
    public List<Interval> splitUpIntervals() {
        return this.intervals;
    }

}
