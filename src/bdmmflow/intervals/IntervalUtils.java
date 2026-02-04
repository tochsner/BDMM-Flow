package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class IntervalUtils {
    public static List<Interval> getIntervals(Parameterization parameterization, double maxIntervalSize) {
        return IntervalUtils.getIntervals(parameterization, maxIntervalSize, 0.0);
    }

    /**
     * Returns a list of intervals that can be used to split up the entire process time. The intervals
     * are built such they don't exceed maxIntervalSize, this means that at least maxProcessTime / maxIntervalSize
     * intervals are returned.
     *
     * @param parameterization - the parameterization object specifying the parameterization intervals.
     * @param maxIntervalSize  - the maximal length of the returned intervals.
     * @return a list of intervals.
     */
    public static List<Interval> getIntervals(Parameterization parameterization, double maxIntervalSize, double rootTime) {
        List<Interval> intervals = new ArrayList<>();

        double currentStartTime = Math.min(0.0, Math.min(rootTime, parameterization.getIntervalEndTimes()[0]));
        double endTime = parameterization.getTotalProcessLength();
        double numIntervals = (int) Math.max(1, Math.ceil((endTime - currentStartTime) / maxIntervalSize));

        for (int i = 0; i < numIntervals; i++) {
            double currentEndTime = currentStartTime + maxIntervalSize;

            // find containing parameterization interval ends
            List<Integer> containingParameterizationIntervalEnds = new ArrayList<>();
            for (int j = 0; j < parameterization.getTotalIntervalCount(); j++) {
                double parameterizationIntervalEndTime = parameterization.getIntervalEndTimes()[j];
                if (
                        Utils.lessThanWithPrecision(currentStartTime, parameterizationIntervalEndTime) &&
                                Utils.lessThanWithPrecision(parameterizationIntervalEndTime, currentEndTime)
                    ) {
                    containingParameterizationIntervalEnds.add(j);
                }
            }


            Interval interval = new Interval(i, containingParameterizationIntervalEnds, currentStartTime, currentEndTime);
            intervals.add(interval);
            currentStartTime += maxIntervalSize;
        }

        return intervals;
    }

    /**
     * Returns the interval corresponding to the given time.
     */
    public static Interval getInterval(double time, List<Interval> intervals) {
        if (time < intervals.get(0).start()) {
            return intervals.get(0);
        }

        for (Interval interval : intervals) {
            if (interval.start() <= time && time <= interval.end()) {
                return interval;
            }
        }

        return intervals.get(intervals.size() - 1);
    }

    public static Set<Integer> getParameterizationIntervalBoundaries(List<Interval> intervals) {
        Set<Integer> boundaries = new HashSet<>();

        return boundaries;
    }
}
