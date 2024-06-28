package bdmmflow.flow.extinctionSystem;

import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 * This class is a lightweight wrapper of the integration output of ExtinctionProbabilitiesODESystem. It allows
 * to conveniently query the extinction probability at a given time.
 */
public class ExtinctionProbabilities {
    protected ContinuousOutputModel[] outputModels;

    public ExtinctionProbabilities(ContinuousOutputModel[] outputModels) {
        this.outputModels = outputModels;
    }

    protected double[] getProbability(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return output.getInterpolatedState();
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
            if (model.getInitialTime() >= time && time >= model.getFinalTime()) {
                return this.getProbability(model, time);
            }
        }

        return this.getProbability(this.outputModels[this.outputModels.length - 1], time);
    }
}
