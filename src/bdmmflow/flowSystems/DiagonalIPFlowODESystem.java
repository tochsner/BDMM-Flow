package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
import org.apache.commons.math3.util.FastMath;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * This class represents the classical flow ODE.
 */
public class DiagonalIPFlowODESystem extends FlowODESystem {

    boolean isNonDiagonalTimeHeterogeneous;

    // cache for matrix exponentials to avoid recomputation
    private final Map<ExponentialCacheKey, RealMatrix> nonDiagonalExpCache = new HashMap<>();

    public DiagonalIPFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance);
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
     * Adds the diagonal time-varying elements to the given system matrix. Note
     * that this only captures the complete behavior if there are no cross-deme
     * birth events.
     */
    void setDiagonalTimeVaryingSystemMatrix(double t, DiagonalMatrix system) {
        double[] extinctProbabilities = this.extinctionProbabilities.getProbability(t);
        int interval = getCurrentParameterizationInterval(t);

        for (int i = 0; i < parameterization.getNTypes(); i++) {
            system.setEntry(
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
            }
        }
    }

    @Override
    protected ContinuousOutputModel integrate(double[] initialState, double start, double end, Interval currentInterval) {
        if (isNonDiagonalTimeHeterogeneous) {
            return super.integrate(initialState, start, end, currentInterval);
        }

        int n = this.parameterization.getNTypes();

        int numSteps = 16;
        double stepSize = (end - start) / (numSteps - 1);

        double[] times = new double[numSteps];
        double[][] states = new double[numSteps][initialState.length];

        // pre-compute non-diagonal exponential (with caching)

        int interval = this.parameterization.getIntervalIndex(start);
        ExponentialCacheKey cacheKey = new ExponentialCacheKey(interval, stepSize);

        RealMatrix nonDiagonalExp = nonDiagonalExpCache.get(cacheKey);
        if (nonDiagonalExp == null) {
            DoubleMatrix constantSystem = Utils.toMatrix(this.buildTimeInvariantSystemMatrix(interval));
            nonDiagonalExp = Utils.toMatrix(MatrixFunctions.expm(
                    constantSystem.mul(stepSize / 2)
            ));
            nonDiagonalExpCache.put(cacheKey, nonDiagonalExp);
        }

        // run Strang scheme

        times[0] = start;
        states[0] = initialState;

        RealMatrix lastStateMatrix = Utils.toMatrix(initialState, n);
        DiagonalMatrix timeVaryingMatrix = new DiagonalMatrix(n);

        for (int i = 1; i < numSteps; i++) {
            double time = start + i*stepSize;
            times[i] = time;

            this.setDiagonalTimeVaryingSystemMatrix(time - 0.5*stepSize, timeVaryingMatrix);

            // compute diagonal exponential

            for (int j = 0; j < n; j++) {
                timeVaryingMatrix.setEntry(j, j, FastMath.exp(timeVaryingMatrix.getEntry(j, j) * stepSize));
            }

            // compute new state using the Strang scheme

            lastStateMatrix = nonDiagonalExp
                    .multiply(timeVaryingMatrix.multiply(nonDiagonalExp)).multiply(lastStateMatrix);
            Utils.fillArray(lastStateMatrix, states[i]);
        }

        return new ContinuousModel(times, states);
    }

    static class ContinuousModel extends ContinuousOutputModel {

        private final double[] times;
        private final double[][] states;

        private double interpolatedTime;

        public ContinuousModel(double[] times, double[][] states) {
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

        public ContinuousModel copy() {
            return new ContinuousModel(this.times, this.states);
        }
    }

    /**
     * cache key for storing matrix exponentials based on interval and step size
     */
    private static class ExponentialCacheKey {
        final int interval;
        final double stepSize;

        ExponentialCacheKey(int interval, double stepSize) {
            this.interval = interval;
            this.stepSize = stepSize;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ExponentialCacheKey that = (ExponentialCacheKey) o;
            return interval == that.interval && Double.compare(that.stepSize, stepSize) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(interval, stepSize);
        }
    }

}
