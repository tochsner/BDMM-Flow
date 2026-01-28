package bdmmflow.intervals;

import java.util.List;

/**
 * This record stores an interval.
 * @param interval - the index of the interval.
 * @param parameterizationIntervals - the indices of the corresponding parameterization intervals.
 * @param start - the start time.
 * @param end - the end time.
 */
public record Interval(int interval, List<Integer> parameterizationIntervals, double start, double end) { }
