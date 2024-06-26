package bdmmflow.flow.intervals;

/**
 * This record stores an interval.
 * @param interval - the index of the interval.
 * @param parameterizationInterval - the index of the corresponding parameterization interval.
 * @param start - the start time.
 * @param end - the end time.
 */
public record Interval(int interval, int parameterizationInterval, double start, double end) { }
