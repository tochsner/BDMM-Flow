package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.utils.LinearTimeInvMatrixSystem;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
import org.apache.commons.math3.util.FastMath;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.List;

/**
 * This class represents the classical flow ODE.
 */
public class DiagonalIPFlowODESystem extends FlowODESystem {

    boolean isNonDiagonalTimeHeterogeneous;

    public DiagonalIPFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, extinctionProbabilities, absoluteTolerance, relativeTolerance);
        isNonDiagonalTimeHeterogeneous = isNonDiagonalTimeHeterogeneous();
    }

    private boolean isNonDiagonalTimeHeterogeneous() {
        for (int i = 0; i < crossBirthRates.length; i++) {
            for (int j = 0; j < crossBirthRates[i].length; j++) {
                for (int k = 0; k < crossBirthRates[i][j].length; k++) {
                    if (crossBirthRates[i][j][k] != 0.0) return true;
                }
            }
        }
        return false;
    }

    /**
     * Integrates over the system backwards in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialStates the initial states at the interval ends.
     * @param intervals the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateBackwards(List<double[]> initialStates, List<Interval> intervals, boolean alwaysStartAtInitialState) {
        if (isNonDiagonalTimeHeterogeneous) {
            return super.integrateBackwards(initialStates, intervals, alwaysStartAtInitialState, false);
        }

        this.currentParameterizationInterval = this.param.getTotalIntervalCount() - 1;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        double[] state = initialStates.get(initialStates.size() - 1).clone();
        for (int i = intervals.size() - 1; i >= 0; i--) {
            Interval interval = intervals.get(i);

            if (alwaysStartAtInitialState) {
                state = initialStates.get(i).clone();
            }

            if (this.currentParameterizationInterval != interval.parameterizationInterval()) {
                assert this.currentParameterizationInterval - 1 == interval.parameterizationInterval();

                this.handleParameterizationIntervalBoundary(
                        interval.end(),
                        this.currentParameterizationInterval,
                        this.currentParameterizationInterval - 1,
                        state
                );
            }

            outputModels[intervals.size() - interval.interval() - 1] = integrate(
                    interval.end(), interval.start(), state
            );
        }

        return outputModels;
    }

    private SimpleContinuousModel integrate(double start, double end, double[] state) {
        int n = this.param.getNTypes();

        int numSteps = 10;
        double stepSize = (end - start) / (numSteps - 1);

        double[] times = new double[numSteps];
        double[][] states = new double[numSteps][state.length];

        // pre-compute non-diagonal exponential

        int interval = this.param.getIntervalIndex(start);
        DoubleMatrix constantSystem = Utils.toMatrix(this.buildTimeInvariantSystemMatrix(interval));
        RealMatrix nonDiagonalExp = Utils.toMatrix(MatrixFunctions.expm(
                constantSystem.mul(stepSize / 2)
        ));

        // run Strang scheme

        times[0] = start;
        states[0] = state;

        RealMatrix lastStateMatrix = Utils.toMatrix(state, n);
        RealMatrix timeVaryingMatrix = new BlockRealMatrix(n, n);

        for (int i = 1; i < numSteps; i++) {
            double time = start + i*stepSize;
            times[i] = time;

            timeVaryingMatrix = timeVaryingMatrix.scalarMultiply(0.0);
            this.addTimeVaryingSystemMatrix(time - 0.5*stepSize, timeVaryingMatrix);

            // compute diagonal exponential

            for (int j = 0; j < n; j++) {
                timeVaryingMatrix.setEntry(j, j, FastMath.exp(timeVaryingMatrix.getEntry(j, j) * stepSize));
            }

            // compute new state using the Strang scheme

            lastStateMatrix = nonDiagonalExp
                    .multiply(timeVaryingMatrix).multiply(nonDiagonalExp).multiply(lastStateMatrix);
            Utils.fillArray(lastStateMatrix, states[i]);
        }

        return new SimpleContinuousModel(times, states);
    }

    private class SimpleContinuousModel extends ContinuousOutputModel {

        private final double[] times;
        private final double[][] states;

        private double interpolatedTime;

        public SimpleContinuousModel(double[] times, double[][] states) {
            this.times = times;
            this.states = states;
        }

        public double getInitialTime() {
            return this.times[0];
        }

        public double getFinalTime() {
            return this.times[this.times.length - 1];
        }

        public double getInterpolatedTime() {
            return interpolatedTime;

        }

        public void setInterpolatedTime(double time) {
            interpolatedTime = time;
        }

        public double[] getInterpolatedState() throws MaxCountExceededException {
            // times are in descending order (backward integration)
            double minTime = times[times.length - 1];
            double maxTime = times[0];

            if (interpolatedTime < minTime - 1e-10 || interpolatedTime > maxTime + 1e-10) {
                throw new OutOfRangeException(interpolatedTime, minTime, maxTime);
            }

            // find the two surrounding time points

            int leftIndex = 0;
            for (int i = 0; i < times.length - 1; i++) {
                if (interpolatedTime <= times[i] && interpolatedTime >= times[i + 1]) {
                    leftIndex = i;
                    break;
                }
            }

            // handle exact match

            if (interpolatedTime == times[leftIndex]) {
                return states[leftIndex].clone();
            }

            int rightIndex = leftIndex + 1;
            if (rightIndex >= times.length) {
                return states[leftIndex].clone();
            }

            // linear interpolation

            double t0 = times[leftIndex];
            double t1 = times[rightIndex];
            double alpha = (interpolatedTime - t0) / (t1 - t0);

            double[] result = new double[states[leftIndex].length];
            for (int i = 0; i < result.length; i++) {
                result[i] = states[leftIndex][i] * (1 - alpha) + states[rightIndex][i] * alpha;
            }

            return result;
        }

        public void append(ContinuousOutputModel model) throws MathIllegalArgumentException, MaxCountExceededException {
            throw new UnsupportedOperationException();
        }

        public void init(double t0, double[] y0, double t) {
            throw new UnsupportedOperationException();
        }

        public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
            throw new UnsupportedOperationException();
        }

        public double[] getInterpolatedDerivatives() throws MaxCountExceededException {
            throw new UnsupportedOperationException();
        }

        public double[] getInterpolatedSecondaryState(int secondaryStateIndex) throws MaxCountExceededException {
            throw new UnsupportedOperationException();
        }

        public double[] getInterpolatedSecondaryDerivatives(int secondaryStateIndex) throws MaxCountExceededException {
            throw new UnsupportedOperationException();
        }
    }

}
