package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class IntervalUtils {
    /**
     * Returns a list of intervals that can be used to split up the entire process time. The intervals
     * are built such they don't exceed maxIntervalSize and that the boundaries of the parameterization
     * intervals are respected.
     * This means that at least maxProcessTime / maxIntervalSize and at least totalParameterizationIntervalCount
     * intervals are returned.
     *
     * @param parameterization - the parameterization object specifying the parameterization intervals.
     * @param maxIntervalSize  - the maximal length of the returned intervals.
     * @return a list of intervals.
     */
    public static List<Interval> getIntervals(Parameterization parameterization, double maxIntervalSize) {
        List<Interval> intervals = new ArrayList<>();

        int currentParameterizationInterval = 0;

        double currentStartTime = 0.0;

        while (currentParameterizationInterval < parameterization.getTotalIntervalCount()) {
            double currentParameterizationIntervalEndTime = parameterization.getIntervalEndTimes()[currentParameterizationInterval];
            double currentMaxEndTime = currentStartTime + maxIntervalSize;

            if (currentParameterizationIntervalEndTime < currentMaxEndTime || Utils.equalWithPrecision(currentParameterizationIntervalEndTime, currentMaxEndTime)) {
                // the current interval goes until the end of the parameterization interval
                intervals.add(
                        new Interval(
                                intervals.size(),
                                currentParameterizationInterval,
                                currentStartTime,
                                currentParameterizationIntervalEndTime
                        )
                );

                currentStartTime = currentParameterizationIntervalEndTime;
                currentParameterizationInterval += 1;
            } else {
                // the current interval is shorter than the current parameterization interval
                intervals.add(
                        new Interval(
                                intervals.size(),
                                currentParameterizationInterval,
                                currentStartTime,
                                currentMaxEndTime
                        )
                );

                currentStartTime = currentMaxEndTime;
            }
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

        for (int i = 1; i < intervals.size(); i++) {
            if (intervals.get(i - 1).parameterizationInterval() != intervals.get(i).parameterizationInterval()) {
                boundaries.add(intervals.get(i).interval());
            }
        }

        return boundaries;
    }
}
