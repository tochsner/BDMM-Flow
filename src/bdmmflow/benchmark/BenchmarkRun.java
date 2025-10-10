package bdmmflow.benchmark;

import java.util.HashMap;
import java.util.Map;

public class BenchmarkRun {
    public static Map<String, String> lastLoggedMetrics = new HashMap<>();

    public static void logMetric(String key, String value) {
        lastLoggedMetrics.put(key, value);
    }

    long duration;
    double likelihood;
    Map<String, String> loggedMetrics;

    public BenchmarkRun(long duration, double likelihood) {
        this.duration = duration;
        this.likelihood = likelihood;

        this.loggedMetrics = this.lastLoggedMetrics;
        this.lastLoggedMetrics = new HashMap<>();
    }
}
