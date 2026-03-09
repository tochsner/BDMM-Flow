package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

public class IPFlow extends Flow {

    Parameterization parameterization;
    final RealMatrix[] timeInvariantQMatrices;
    final RealMatrix[] timeInvariantQTMatrices;
    final RealMatrix[] timeInvariantTMatrices;

    public IPFlow(ContinuousOutputModel[] outputModels, int n, List<InitialState> initialStateArrays, boolean wasInitialStateResetAtEachInterval, Parameterization parameterization, RealMatrix[] timeInvariantQMatrices, RealMatrix[] timeInvariantQTMatrices, RealMatrix[] timeInvariantTMatrices) {
        super(outputModels, n, initialStateArrays, wasInitialStateResetAtEachInterval);
        this.parameterization = parameterization;
        this.timeInvariantQMatrices = timeInvariantQMatrices;
        this.timeInvariantQTMatrices = timeInvariantQTMatrices;
        this.timeInvariantTMatrices = timeInvariantTMatrices;
    }

    RealMatrix getFlow(int interval, double time) {
        ContinuousOutputModel output = this.outputModels[interval];

        double[] state;
        synchronized (output) {
            output.setInterpolatedTime(time);
            state = output.getInterpolatedState().clone();
        }
        RealMatrix flow = Utils.toMatrix(state, this.n);

        double negDeltaT = time - this.outputModels[interval].getInitialTime();
        RealMatrix Q = this.timeInvariantQMatrices[interval];
        RealMatrix QT = this.timeInvariantQTMatrices[interval];
        RealMatrix T = this.timeInvariantTMatrices[interval];

        RealMatrix expT = QT.multiply(Utils.expmUpperTriangular(T.scalarMultiply(negDeltaT))).multiply(Q);

        synchronized (this) {
            a++;
            b += new SingularValueDecomposition(flow).getConditionNumber();
        }

        return expT.multiply(flow);
    }

}
