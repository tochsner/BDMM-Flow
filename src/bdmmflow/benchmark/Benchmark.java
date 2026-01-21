package bdmmflow.benchmark;

import bdmmflow.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdmmprime.trajectories.simulation.SimulatedTree;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;

import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

public class Benchmark {

    public static void main(String[] args) {
        int NUM_TRIALS = 30_000;

        ParameterizationSampler sampler = new ParameterizationSampler();

        List<BenchmarkResult> results = runBenchmarks(NUM_TRIALS, sampler);
        writeResults(results, "results.csv");
    }

    static List<BenchmarkResult> runBenchmarks(int numTrials, ParameterizationSampler sampler) {
        List<BenchmarkResult> results = new LinkedList<>();

        for (int i = 0; i < numTrials; i++) {
            Parameterization parameterization = sampler.sampleParameterization();
            RealParameter startTypePriorProbs = sampler.sampleStartTypePriorProbs(parameterization);
            Tree tree = simulateTree(parameterization, startTypePriorProbs);

            BenchmarkRun bdmmRun = runBDMMBenchmark(tree, parameterization, startTypePriorProbs);

            String[] initialStateStrategies = new String[]{
                    "identity",
                    "random",
            };

            for (String strategy : initialStateStrategies) {
                boolean useSplitting = false;
                boolean useInverseFlow = true;
                BenchmarkRun flowRun = runFlowBenchmark(tree, parameterization, startTypePriorProbs, useInverseFlow, useSplitting, strategy);
                BenchmarkResult result = new BenchmarkResult(
                        i, parameterization, tree, flowRun, bdmmRun, useInverseFlow, useSplitting, strategy
                );
                results.add(result);

                useInverseFlow = false;
                flowRun = runFlowBenchmark(tree, parameterization, startTypePriorProbs, useInverseFlow, useSplitting, strategy);
                result = new BenchmarkResult(
                        i, parameterization, tree, flowRun, bdmmRun, useInverseFlow, useSplitting, strategy
                );
                results.add(result);

                useSplitting = true;
                flowRun = runFlowBenchmark(tree, parameterization, startTypePriorProbs, useInverseFlow, useSplitting, strategy);
                result = new BenchmarkResult(
                        i, parameterization, tree, flowRun, bdmmRun, useInverseFlow, useSplitting, strategy
                );
                results.add(result);
            }

            if (i % 100 == 0) {
                System.out.println(i);
                writeResults(results, "results.csv");
            }
        }

        return results;
    }

    static Tree simulateTree(Parameterization parameterization, RealParameter startTypePriorProbs) throws IllegalStateException {
        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", new RealParameter("0.0"),
                "startTypePriorProbs", startTypePriorProbs,
                "minSamples", 2,
                "simulateUntypedTree", true
        );
        return simulatedTree;
    }

    static BenchmarkRun runFlowBenchmark(
            Tree tree,
            Parameterization parameterization,
            RealParameter startTypePriorProbs,
            boolean useInverseFlow,
            boolean useSplitting,
            String initialStateStrategy
    ) {
        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "startTypePriorProbs", startTypePriorProbs,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "useODESplitting", useSplitting,
                "parallelize", true
        );
        density.initAndValidate();

        long start = System.nanoTime();
        double likelihood = density.calculateLogP();
        long duration = System.nanoTime() - start;

        return new BenchmarkRun(duration, likelihood);
    }

    static BenchmarkRun runBDMMBenchmark(Tree tree, Parameterization parameterization, RealParameter startTypePriorProbs) {
        bdmmprime.distribution.BirthDeathMigrationDistribution density = new bdmmprime.distribution.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "startTypePriorProbs", startTypePriorProbs,
                "typeLabel", "type",
                "parallelize", true
        );
        density.initAndValidate();

        long start = System.nanoTime();
        double likelihood = density.calculateLogP();
        long duration = System.nanoTime() - start;

        return new BenchmarkRun(duration, likelihood);
    }

    static void writeResults(List<BenchmarkResult> results, String fileName) {
        try (FileWriter fileWriter = new FileWriter(fileName)) {
            fileWriter.write(results.get(0).getHeaders());
            fileWriter.write("\n");

            for (BenchmarkResult result : results) {
                fileWriter.write(result.toString());
                fileWriter.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        LoggedMetric.storeMetrics("metrics.csv");
    }

}
