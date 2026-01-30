package bdmmflow.extinctionSystem;

import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 * This class is a lightweight wrapper of the integration output of ExtinctionProbabilitiesODESystem. It allows
 * to conveniently query the extinction probability at a given time.
 */
public class ExtinctionProbabilities {
    ContinuousOutputModel[] outputModels;

    public ExtinctionProbabilities(ContinuousOutputModel[] outputModels) {
        this.outputModels = outputModels;
    }

    public double[] getProbability(ContinuousOutputModel output, double time) {
        synchronized (output) {
            output.setInterpolatedTime(time);
            return output.getInterpolatedState().clone();
        }
    }

    public double[] unsafeGetProbability(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return output.getInterpolatedState();
    }

    /**
     * Returns the extinction probability at the given time.
     * @param time the time.
     * @return the extinction probability at the given time.
     */
    public double[] getProbability(double time) {
        return this.getProbability(this.getOutputModel(time), time);
    }

    /**
     * Returns the continuous output model for the given time.
     */
    public ContinuousOutputModel getOutputModel(double time) {
        if (time > this.outputModels[0].getInitialTime()) {
            return this.outputModels[0];
        }

        for (ContinuousOutputModel model : this.outputModels) {
            if (model.getInitialTime() >= time && time > model.getFinalTime()) {
                return model;
            }
        }

        return this.outputModels[this.outputModels.length - 1];
    }

}
