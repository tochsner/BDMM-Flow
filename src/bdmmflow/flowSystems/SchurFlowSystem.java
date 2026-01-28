package bdmmflow.flowSystems;
import bdmmflow.utils.Utils;
import com.flag4j.Matrix;
import com.flag4j.linalg.decompositions.RealSchurDecomposition;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

public class SchurFlowSystem extends FlowODESystem {

    RealMatrix A_0, U, U_T, T;

    public SchurFlowSystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance, int seed, double maxConditionNumber) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance, seed, maxConditionNumber);

        double endTime = intervals.get(intervals.size() - 1).end();
        this.A_0 = this.buildSystemMatrix(endTime);
        Matrix initialSystemMatrixFlag4J = Utils.toFlag4JMatrix(this.A_0);

        RealSchurDecomposition realSchurDecomposition = new RealSchurDecomposition();
        realSchurDecomposition.decompose(initialSystemMatrixFlag4J);

        this.U = Utils.toMatrix(realSchurDecomposition.getRealU());
        this.U_T = this.U.transpose();
        this.T = Utils.toMatrix(realSchurDecomposition.getRealT());
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int n = this.parameterization.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, n);
        RealMatrix currentA = this.buildSystemMatrix(t);

        RealMatrix systemMatrix = this.T.add(this.U_T.multiply(currentA.subtract(this.A_0)).multiply(this.U));

        RealMatrix yDotMatrix = systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
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

        return new SchurFlow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalsBoundaries,
                this.U
        );
    }
}
