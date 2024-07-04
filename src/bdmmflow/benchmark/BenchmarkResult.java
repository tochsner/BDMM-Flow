package bdmmflow.benchmark;

import bdmmprime.parameterization.Parameterization;
import beast.base.evolution.tree.Tree;

import java.util.StringJoiner;


public class BenchmarkResult {

    Parameterization parameterization;
    Tree tree;
    BenchmarkRun flowRun;
    BenchmarkRun bdmmRun;
    boolean useInverseFlow;
    boolean useRandomInitialMatrix;
    int minNumIntervals;

    public BenchmarkResult(
            Parameterization parameterization,
            Tree tree,
            BenchmarkRun flowRun,
            BenchmarkRun bdmmRun,
            boolean useInverseFlow,
            boolean useRandomInitialMatrix,
            int minNumIntervals
    ) {
        this.parameterization = parameterization;
        this.tree = tree;
        this.flowRun = flowRun;
        this.bdmmRun = bdmmRun;
        this.useInverseFlow = useInverseFlow;
        this.useRandomInitialMatrix = useRandomInitialMatrix;
        this.minNumIntervals = minNumIntervals;
    }

    @Override
    public String toString() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add(Integer.toString(this.tree.getNodeCount()));
        joiner.add(Integer.toString(this.tree.getLeafNodeCount()));

        joiner.add(Integer.toString(this.parameterization.getNTypes()));
        joiner.add(Double.toString(this.parameterization.getTotalProcessLength()));

        joiner.add(Double.toString(this.flowRun.likelihood));
        joiner.add(Long.toString(this.flowRun.duration));

        joiner.add(Double.toString(this.bdmmRun.likelihood));
        joiner.add(Long.toString(this.bdmmRun.duration));

        joiner.add(Boolean.toString(this.useInverseFlow));
        joiner.add(Boolean.toString(this.useRandomInitialMatrix));
        joiner.add(Integer.toString(this.minNumIntervals));

        return joiner.toString();
    }

    public static String getHeaders() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add("node_count");
        joiner.add("leaf_count");

        joiner.add("types_count");
        joiner.add("process_length");

        joiner.add("flow_likelihood");
        joiner.add("flow_duration");

        joiner.add("bdmm_likelihood");
        joiner.add("bdmm_duration");

        joiner.add("use_inverse_flow");
        joiner.add("use_random_initial_matrix");
        joiner.add("num_intervals");

        return joiner.toString();
    }
}
