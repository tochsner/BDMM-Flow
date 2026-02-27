package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

public class SchurFlow extends Flow {

    RealMatrix U;

    public SchurFlow(ContinuousOutputModel[] outputModels, int n, List<double[]> initialStateArrays, boolean wasInitialStateResetAtEachInterval, RealMatrix U) {
        super(outputModels, n, initialStateArrays, wasInitialStateResetAtEachInterval);
        this.U = U;
    }

    RealMatrix getFlow(ContinuousOutputModel output, double time) {
        synchronized (output) {
            output.setInterpolatedTime(time);
            return this.U.multiply(Utils.toMatrix(output.getInterpolatedState(), this.n));
        }
    }

}
