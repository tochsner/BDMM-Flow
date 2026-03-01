package bdmmflow.flowSystems;

import org.apache.commons.math3.linear.RealMatrix;

public record InitialState(
        double[] initialState,
        RealMatrix inverse
) { }
