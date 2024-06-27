package bdmmflow.flow.flowSystems;

import bdmmflow.flow.Utils;
import bdmmflow.flow.extinctionSystem.ExtinctionProbabilities;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class Flow extends ExtinctionProbabilities {

    RealMatrix inverseInitialState;
    boolean wasInitialStateResetAtEachInterval;

    int n;

    public Flow(ContinuousOutputModel[] flows, int n) {
        this(flows, n, MatrixUtils.createRealIdentityMatrix(n), false);
    }

    public Flow(ContinuousOutputModel[] flows, int n, RealMatrix inverseInitialState, boolean wasInitialStateResetAtEachInterval) {
        super(flows);

        this.n = n;
        this.inverseInitialState = inverseInitialState;
        this.wasInitialStateResetAtEachInterval = wasInitialStateResetAtEachInterval;
    }

    protected RealMatrix getFlow(ContinuousOutputModel output, double time) {
        return Utils.toMatrix(super.getProbability(output, time), this.n);
    }

    public int getInterval(double time) {
        for (int i = 0; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                return i;
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }

    public RealMatrix getFlow(double time, int startingAtInterval) {
        RealMatrix accumulatedFlow = null;

        for (int i = startingAtInterval; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                if (accumulatedFlow == null) {
                    return this.getFlow(this.outputModels[i], time);
                } else {
                    return this.getFlow(this.outputModels[i], time).multiply(accumulatedFlow);
                }
            }

            if (this.wasInitialStateResetAtEachInterval) {
                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
                accumulatedFlow = this.inverseInitialState.multiply(flowEnd.multiply(accumulatedFlow));
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }
}