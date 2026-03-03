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

    public static List<Interval> getIntervals(Parameterization parameterization, double maxIntervalSize, double rootTime) {
        List<Interval> intervals = new ArrayList<>();

        int currentParameterizationInterval = 0;

        double currentStartTime = Math.min(0.0, Math.min(rootTime, parameterization.getIntervalEndTimes()[0]));

        while (currentParameterizationInterval < parameterization.getTotalIntervalCount()) {
            double currentParameterizationIntervalEndTime = parameterization.getIntervalEndTimes()[currentParameterizationInterval];

            if (Utils.equalWithPrecision(currentParameterizationIntervalEndTime, currentStartTime) || currentParameterizationIntervalEndTime < currentStartTime) {
                // the current interval is empty or ends before it starts. this can happen when the interval has a negative end time.
                currentParameterizationInterval += 1;
                continue;
            }

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
