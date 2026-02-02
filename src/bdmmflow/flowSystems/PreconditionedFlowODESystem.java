package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalUtils;
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
public class PreconditionedFlowODESystem extends FlowODESystem {
    RealMatrix[] preconditioners;
    RealMatrix[] inversePreconditioners;

    public PreconditionedFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance, int seed, double maxConditionNumber) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance, seed, maxConditionNumber);
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int numTypes = this.parameterization.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);

        Interval interval = IntervalUtils.getInterval(t, this.intervals);
        RealMatrix preconditioner = this.preconditioners[interval.interval()];
        RealMatrix inversePreconditioner = this.inversePreconditioners[interval.interval()];

        RealMatrix yDotMatrix = inversePreconditioner.multiply(systemMatrix).multiply(preconditioner).multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }

    List<double[]> getInitialStates(String initialMatrixStrategy, List<Interval> intervals) {
        return switch (initialMatrixStrategy) {
            case "random" -> {
                RealMatrix matrix = Utils.getRandomMatrix(this.parameterization.getNTypes(), this.seed);

                List<double[]> arrays = new ArrayList<>();
                for (Interval interval : intervals) {
                    RealMatrix inversePreconditioner = this.inversePreconditioners[interval.interval()];

                    double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                    Utils.fillArray(inversePreconditioner.multiply(matrix), array);

                    arrays.add(array);
                }

                yield arrays;
            }
            case "identity" -> {
                RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());

                List<double[]> arrays = new ArrayList<>();
                for (Interval interval : intervals) {
                    RealMatrix inversePreconditioner = this.inversePreconditioners[interval.interval()];

                    double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                    Utils.fillArray(inversePreconditioner.multiply(matrix), array);

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
        this.computePreconditioners(intervals);
        List<double[]> initialStates = this.getInitialStates(initialMatrixStrategy, intervals);

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStates,
                intervals,
                resetInitialStateAtIntervalsBoundaries,
                parallelize
        );

        return new PreconditionedFlow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries,
                preconditioners,
                inversePreconditioners,
                intervals
        );
    }

    private void computePreconditioners(List<Interval> intervals) {
        this.preconditioners = new RealMatrix[intervals.size()];
        this.inversePreconditioners = new RealMatrix[intervals.size()];

        int n = this.parameterization.getNTypes();

        for (int i = 0; i < intervals.size(); i++) {
            RealMatrix preconditioner = new BlockRealMatrix(n, n);
            RealMatrix systemMatrix = this.buildSystemMatrix(intervals.get(i).end());

            preconditioner = systemMatrix;

//            for (int j = 0; j < n; j++) {
//                for (int k = 0; k < n; k++) {
//                    if (j != k) continue;
//                    preconditioner.setEntry(j, k, systemMatrix.getEntry(j, k));
//                }
//            }

            this.preconditioners[intervals.get(i).interval()] = preconditioner;
            this.inversePreconditioners[intervals.get(i).interval()] = MatrixUtils.inverse(preconditioner);
        }
    }
}
