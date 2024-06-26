package bdmmflow.flow;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class ContinuousIntervalOutputModel {
    ContinuousOutputModel[] outputModels;

    public ContinuousIntervalOutputModel(ContinuousOutputModel[] flows) {
        this.outputModels = flows;
    }

    protected double[] getOutput(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return output.getInterpolatedState();
    }

    public double[] getOutput(double time) {
        for (ContinuousOutputModel model : this.outputModels) {
            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                return this.getOutput(model, time);
            }
        }

        if (time > this.outputModels[0].getInitialTime()) {
            return this.getOutput(this.outputModels[0], time);
        } else {
            return this.getOutput(this.outputModels[this.outputModels.length - 1], time);
        }
    }
}
