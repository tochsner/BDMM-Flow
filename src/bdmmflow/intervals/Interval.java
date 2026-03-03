package bdmmflow.intervals;

/**
 * This record stores an interval.
 * @param interval - the index of the interval.
 * @param start - the start time.
 * @param end - the end time.
 */
public record Interval(int interval, double start, double end) { }
