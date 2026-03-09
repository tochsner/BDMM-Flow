package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.hipparchus.linear.SchurTransformer;

import java.util.List;

/**
 * This class represents the classical flow ODE.
 */
public class IPFlowODESystem extends FlowODESystem {

    final RealMatrix[] timeInvariantQMatrices;
    final RealMatrix[] timeInvariantQTMatrices;
    final RealMatrix[] timeInvariantTMatrices;

    public IPFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities, List<Interval> intervals, double absoluteTolerance, double relativeTolerance, int seed, double maxConditionNumber, boolean useLoucaPennellIntervals) {
        super(parameterization, extinctionProbabilities, intervals, absoluteTolerance, relativeTolerance, seed, maxConditionNumber, useLoucaPennellIntervals);

        this.timeInvariantQMatrices = new RealMatrix[intervals.size()];
        this.timeInvariantQTMatrices = new RealMatrix[intervals.size()];
        this.timeInvariantTMatrices = new RealMatrix[intervals.size()];
        for (int i = 0; i < intervals.size(); i++) {
            RealMatrix invariantMatrix = this.timeInvariantSystemMatrices[i];
            SchurTransformer schur = new SchurTransformer(Utils.toHipparchusMatrix(invariantMatrix));
            this.timeInvariantQMatrices[i] = Utils.toMatrix(schur.getP());
            this.timeInvariantQTMatrices[i] = this.timeInvariantQMatrices[i].transpose();
            this.timeInvariantTMatrices[i] = Utils.toMatrix(schur.getT());
        }
    }

    /**
     * Builds the system matrix for a given time point.
     */
    RealMatrix buildSystemMatrix(double t) {
//        RealMatrix timeInvariantMatrix = this.buildTimeInvariantSystemMatrix(t);
//        RealMatrix timeVaryingMatrix = this.buildTimeVaryingSystemMatrix(t);
//
//        double[] values = new double[this.param.getNTypes() * this.param.getNTypes()];
//        Utils.fillArray(timeVaryingMatrix, values);
//
//        RealMatrix A = Utils.toMatrix(MatrixFunctions.expm(this.jblasConstantMatrix.mul(this.param.getTotalProcessLength() - t)));
//        RealMatrix B = Utils.toMatrix(MatrixFunctions.expm(this.jblasConstantMatrix.mul(t - this.param.getTotalProcessLength())));
//
//        return A.multiply(timeVaryingMatrix.add(timeInvariantMatrix).subtract(constantMatrix)).multiply(B);

        int interval = this.parameterization.getIntervalIndex(t);
        double deltaT = this.intervals.get(interval).end() - t;

        RealMatrix timeVaryingMatrix = new BlockRealMatrix(this.parameterization.getNTypes(), this.parameterization.getNTypes());
        this.addTimeVaryingSystemMatrix(t, timeVaryingMatrix);

        RealMatrix Q = this.timeInvariantQMatrices[interval];
        RealMatrix QT = this.timeInvariantQTMatrices[interval];
        RealMatrix T = this.timeInvariantTMatrices[interval];

        RealMatrix expT = Q.multiply(Utils.expmUpperTriangular(T.scalarMultiply(deltaT))).multiply(QT);
        RealMatrix expNegT = QT.multiply(Utils.expmUpperTriangular(T.scalarMultiply(-deltaT))).multiply(Q);

        return expT.multiply(timeVaryingMatrix).multiply(expNegT);
    }

    /**
     * Calculates the flow integral using the given intervals.
     *
     * @return the calculated flow.
     */
    @Override
    public IFlow calculateFlowIntegral(
            String initialMatrixStrategy,
            boolean parallelize
    ) {
        // this.splitUpIntervals();
        List<InitialState> initialStates = this.getInitialStates(initialMatrixStrategy, intervals);

        boolean resetInitialStateAtIntervalBoundaries = 1 < this.intervals.size();

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStates.stream().map(InitialState::initialState).toList(),
                intervals,
                resetInitialStateAtIntervalBoundaries,
                parallelize
        );

        return new IPFlow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalBoundaries,
                this.parameterization,
                this.timeInvariantQMatrices,
                this.timeInvariantQTMatrices,
                this.timeInvariantTMatrices
        );
    }

    @Override
    public int getDimension() {
        return parameterization.getNTypes() * parameterization.getNTypes();
    }

}
