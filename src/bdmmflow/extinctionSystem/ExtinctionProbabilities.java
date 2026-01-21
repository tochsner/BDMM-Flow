package bdmmflow.extinctionSystem;

import bdmmflow.intervals.Interval;
import bdmmprime.util.Utils;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.ArrayList;
import java.util.List;

/**
 * This class is a lightweight wrapper of the integration output of ExtinctionProbabilitiesODESystem. It allows
 * to conveniently query the extinction probability at a given time.
 */
public class ExtinctionProbabilities {
    ContinuousOutputModel[] outputModels;

    public ExtinctionProbabilities(ContinuousOutputModel[] outputModels) {
        this.outputModels = outputModels;
    }

    double[] getProbability(ContinuousOutputModel output, double time) {
        synchronized (output) {
            output.setInterpolatedTime(time);
            return output.getInterpolatedState();
        }
    }

    /**
     * Returns the extinction probability at the given time.
     * @param time the time.
     * @return the extinction probability at the given time.
     */
    public double[] getProbability(double time) {
        if (time > this.outputModels[0].getInitialTime()) {
            return this.getProbability(this.outputModels[0], time);
        }

        for (ContinuousOutputModel model : this.outputModels) {
            if (model.getInitialTime() >= time && time > model.getFinalTime()) {
                return this.getProbability(model, time);
            }
        }

        return this.getProbability(this.outputModels[this.outputModels.length - 1], time);
    }

    public ExtinctionProbabilities copy() {
        ContinuousOutputModel[] clonedOutputModels = new ContinuousOutputModel[this.outputModels.length];
        System.arraycopy(this.outputModels, 0, clonedOutputModels, 0, this.outputModels.length);

        return new ExtinctionProbabilities(clonedOutputModels);
    }

}
