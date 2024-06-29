package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;

import java.util.ArrayList;
import java.util.List;

public class IntervalUtils {
    /**
     * Returns a list of intervals that can be used to split up the entire process time. The intervals
     * are built such they don't exceed maxIntervalSize and that the boundaries of the parameterization
     * intervals are respected.
     * This means that at least maxProcessTime / maxIntervalSize and at least totalParameterizationIntervalCount
     * intervals are returned.
     *
     * @param parameterization - the parameterization object specifying the parameterization intervals.
     * @param maxIntervalSize - the maximal length of the returned intervals.
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
}
