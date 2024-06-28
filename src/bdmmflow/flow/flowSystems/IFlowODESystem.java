package bdmmflow.flow.flowSystems;

import bdmmflow.flow.intervals.Interval;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;

public interface IFlowODESystem {
    IFlow calculateFlowIntegral(
            List<Interval> intervals,
            RealMatrix initialState,
            RealMatrix inverseInitialState,
            boolean resetInitialStateAtIntervalsBoundaries
    );
}
