package bdmmflow.flowSystems;

import org.apache.commons.math3.linear.RealMatrix;

public interface IFlow {
    int getInterval(double time);
    RealMatrix getFlow(double time, int startingAtInterval);
    double[] integrateUsingFlow(
            double timeStart,
            double timeEnd,
            double[] endState
    );
}
