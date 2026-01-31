package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.util.FastMath;

import java.util.List;

/**
 * This class is a lightweight wrapper of the result of the Flow ODE integration. It allows to easily query the
 * flow at every point in time and to use the flow to efficiently integrate over a time span.
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class QRFlow extends Flow {

    public QRFlow(ContinuousOutputModel[] outputModels, int n, List<double[]> initialStateArrays, boolean wasInitialStateResetAtEachInterval) {
        super(outputModels, n, initialStateArrays, wasInitialStateResetAtEachInterval);
    }

    /**
     * Allows to integrate over an edge of a tree using the pre-computed flow.
     *
     * @param timeStart the time of the node closer to the root.
     * @param timeEnd   the time of the node closer to the leaves.
     * @param endState  the initial state at the node closer to the leaves.
     * @return the integration result at the time of the node closer to the root.
     */
    @Override
    public double[] integrateUsingFlow(double timeStart, double timeEnd, double[] endState) {
        RealMatrix flowMatrixStart = this.getFlow(timeStart);
        double[] ODEStateEnd = this.getState(timeEnd);

        RealVector likelihoodVectorEnd = Utils.toVector(endState);
        RealVector solution = this.solve(ODEStateEnd, likelihoodVectorEnd, this.n);
        RealVector likelihoodVectorStart = flowMatrixStart.operate(solution);

        return likelihoodVectorStart.toArray();
    }

    RealMatrix getFlow(double time) {
        double[] state = this.getState(time);

        RealMatrix Q = Utils.toMatrix(state, this.n);
        RealMatrix R = Utils.toMatrix(state, this.n, this.n*this.n);

        return Q.multiply(R);
    }

    public double[] getState(double time) {
        int timeInterval = this.getInterval(time);
        return this.getState(this.outputModels[timeInterval], time);
    }

    double[] getState(ContinuousOutputModel output, double time) {
        synchronized (output) {
            output.setInterpolatedTime(time);
            return output.getInterpolatedState();
        }
    }

    public RealVector solve(double[] yState, RealVector b, int n) {
        // Step 1: Compute y = Q^T * b
        // Q is stored column-major: Q[row,col] = yState[col*n + row]
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                // y[i] = (Q^T * b)[i] = sum_j Q[j,i] * b[j]
                sum += yState[i * n + j] * b.getEntry(j);
            }
            y[i] = sum;
        }

        // Step 2: Back-substitution for Rx = y
        // R is stored as full n*n matrix in column-major format starting at yState[n*n]
        // R[row,col] = yState[n*n + col*n + row]
        double[] x = new double[n];
        int rOffset = n * n;

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                // R[i,j] in column-major: yState[rOffset + j*n + i]
                sum += yState[rOffset + j * n + i] * x[j];
            }

            // R[i,i] in column-major: yState[rOffset + i*n + i]
            double rDiag = yState[rOffset + i * n + i];

            x[i] = (y[i] - sum) / rDiag;
        }

        return new ArrayRealVector(x, false);
    }

}
