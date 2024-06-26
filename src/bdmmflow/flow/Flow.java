package bdmmprime.flow;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.util.Pair;

import java.util.LinkedHashMap;
import java.util.Map;

public class Flow extends ContinuousIntervalOutputModel {

    class LRUCache<K, V> extends LinkedHashMap<K, V> {
        private final int cacheSize;

        public LRUCache(int cacheSize) {
            super(20, 0.75F, true);
            this.cacheSize = cacheSize;
        }

        protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
            return size() >= cacheSize;
        }
    }

    LRUCache<Pair<Double, Integer>, RealMatrix> cache = new LRUCache<>(20);

    RealMatrix[] initialFlows;
    RealMatrix inverseInitialState;
    boolean useIntervals;

    int n;
    double[] logScaleFactors;

    public Flow(ContinuousOutputModel[] flows, int n) {
        this(flows, n, MatrixUtils.createRealIdentityMatrix(n), false);
    }

    public Flow(ContinuousOutputModel[] flows, int n, RealMatrix inverseInitialState, boolean useIntervals) {
        super(flows);

        this.n = n;
        this.inverseInitialState = inverseInitialState;
        this.useIntervals = useIntervals;
    }

    protected RealMatrix getFlow(ContinuousOutputModel output, double time) {
        return Utils.toMatrix(super.getOutput(output, time), this.n);
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

    public RealMatrix getFlow(double time, int startingAtInterval, boolean checkCache, boolean storeInCache) {
        Pair<Double, Integer> key = new Pair<>(time, startingAtInterval);

        if (checkCache && this.cache.containsKey(key))
            return this.cache.get(key);

        RealMatrix flow = MatrixUtils.createRealIdentityMatrix(this.n);

        for (int i = startingAtInterval; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (
                    model.getInitialTime() <= time && time <= model.getFinalTime() || model.getFinalTime() <= time && time <= model.getInitialTime()
            ) {
                RealMatrix result = i > startingAtInterval ? flow.multiply(this.getFlow(this.outputModels[i], time)) : this.getFlow(this.outputModels[i], time);
                if (storeInCache) this.cache.put(key, result);

                return result;
            }

            if (this.useIntervals) {
                RealMatrix flowEnd = this.getFlow(this.outputModels[i], this.outputModels[i].getFinalTime());
                flow = flow.multiply(flowEnd).multiply(this.inverseInitialState);
            }
        }

        throw new IllegalArgumentException("The provided time is out of bounds for the stored flows.");
    }
}
