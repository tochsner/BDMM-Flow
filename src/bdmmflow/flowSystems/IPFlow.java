package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

public class IPFlow extends Flow {

    Parameterization parameterization;

    public IPFlow(ContinuousOutputModel[] outputModels, int n, List<InitialState> initialStateArrays, boolean wasInitialStateResetAtEachInterval, Parameterization parameterization) {
        super(outputModels, n, initialStateArrays, wasInitialStateResetAtEachInterval);
        this.parameterization = parameterization;
    }

    RealMatrix getFlow(int interval, double time) {
        ContinuousOutputModel output = this.outputModels[interval];

        double[] state;
        synchronized (output) {
            output.setInterpolatedTime(time);
            state = output.getInterpolatedState().clone();
        }
        RealMatrix flow = Utils.toMatrix(state, this.n);

        double c = new SingularValueDecomposition(flow).getConditionNumber();

        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                double factor = Math.exp(state[this.n*this.n + i]);
                flow.multiplyEntry(i, j, factor);
            }
        }

        double d = new SingularValueDecomposition(flow).getConditionNumber();

        synchronized (this) {
            if (Double.isFinite(c)) {
                a++;
                b += c / d;
                System.out.println(c / d + "(" + c + ") running " + b / a);
            }
        }



        return flow;
    }

}
