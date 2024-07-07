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
import java.util.Random;

public class Benchmark {

    public static void main(String[] args) {
        int NUM_TIMES = 30000;

        ParameterizationSampler sampler = new ParameterizationSampler();

        List<BenchmarkResult> results = runBenchmarks(NUM_TIMES, sampler);
        writeResults(results, "results.csv");
    }

    static List<BenchmarkResult> runBenchmarks(int numTrials, ParameterizationSampler sampler) {
        List<BenchmarkResult> results = new LinkedList<>();

        for (int i = 0; i < numTrials; i++) {
            try {
                Parameterization parameterization = sampler.sampleParameterization();
                RealParameter frequencies = sampler.sampleFrequencies(parameterization);
                Tree tree = simulateTree(parameterization, frequencies);

                boolean useInverseFlow = Math.random() > 0.5;
                boolean useRandom = Math.random() > 0.5;
                int minNumIntervals = Math.random() > 0.5 ? 1 : 4;

                BenchmarkRun flowRun = runFlowBenchmark(tree, parameterization, frequencies, useInverseFlow, useRandom, minNumIntervals);
                BenchmarkRun bdmmRun = runBDMMBenchmark(tree, parameterization, frequencies);

                BenchmarkResult result = new BenchmarkResult(parameterization, tree, flowRun, bdmmRun, useInverseFlow, useRandom, minNumIntervals);
                results.add(result);
            } catch (RuntimeException ignored) {
                System.out.println(ignored);
            }

            if (i % 100 == 0) {
                System.out.println(i);
                writeResults(results, "results.csv");
            }
        }

        return results;
    }

    static Tree simulateTree(Parameterization parameterization, RealParameter frequencies) throws IllegalStateException {
        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", new RealParameter("0.0"),
                "frequencies", frequencies,
                "minSamples", 2,
                "simulateUntypedTree", true
        );
        return simulatedTree;
    }

    static BenchmarkRun runFlowBenchmark(
            Tree tree,
            Parameterization parameterization,
            RealParameter frequencies,
            boolean useInverseFlow,
            boolean useRandomInitialMatrix,
            int minNumIntervals
    ) {
        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "frequencies", frequencies,
                "typeLabel", "type",
                "useRandomInitialMatrix", useRandomInitialMatrix,
                "useInverseFlow", useInverseFlow,
                "minNumIntervals", minNumIntervals
        );
        density.initAndValidate();

        long start = System.nanoTime();
        double likelihood = density.calculateLogP();
        long duration = System.nanoTime() - start;

        return new BenchmarkRun(duration, likelihood);
    }

    static BenchmarkRun runBDMMBenchmark(Tree tree, Parameterization parameterization, RealParameter frequencies) {
        bdmmprime.distribution.BirthDeathMigrationDistribution density = new bdmmprime.distribution.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "frequencies", frequencies,
                "typeLabel", "type",
                "parallelize", false
        );
        density.initAndValidate();

        long start = System.nanoTime();
        double likelihood = density.calculateLogP();
        long duration = System.nanoTime() - start;

        return new BenchmarkRun(duration, likelihood);
    }

    static void writeResults(List<BenchmarkResult> results, String fileName) {
        try (FileWriter fileWriter = new FileWriter(fileName)) {
            fileWriter.write(BenchmarkResult.getHeaders());
            fileWriter.write("\n");

            for (BenchmarkResult result : results) {
                fileWriter.write(result.toString());
                fileWriter.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}
