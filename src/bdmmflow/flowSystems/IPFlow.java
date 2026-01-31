package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.jblas.MatrixFunctions;

import java.util.List;

public class IPFlow extends Flow {

    Parameterization parameterization;

    public IPFlow(ContinuousOutputModel[] outputModels, int n, List<double[]> initialStateArrays, boolean wasInitialStateResetAtEachInterval, Parameterization parameterization) {
        super(outputModels, n, initialStateArrays, wasInitialStateResetAtEachInterval);
        this.parameterization = parameterization;
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        double[] state;
        synchronized (output) {
            output.setInterpolatedTime(time);
            state = output.getInterpolatedState();
        }
        RealMatrix flow = Utils.toMatrix(state, this.n);

        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                double factor = Math.exp(state[this.n*this.n + i]);
                flow.multiplyEntry(i, j, factor);
            }
        }

        return flow;
    }

}
