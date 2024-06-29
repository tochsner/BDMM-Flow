package bdmmflow.utils;

import org.apache.commons.math3.linear.*;

import java.util.Arrays;
import java.util.Random;

public class Utils {

    private static final Random random = new Random();

    public static boolean isSingular(RealMatrix matrix) {
        if (matrix.getRowDimension() != matrix.getColumnDimension()) return false;

        SingularValueDecomposition svd = new SingularValueDecomposition(matrix);
        return svd.getRank() != matrix.getColumnDimension();
    }

    public static RealMatrix getRandomMatrix(int dimension, int seed) {
        RealMatrix randomMatrix;

        do {
            randomMatrix = new BlockRealMatrix(dimension, dimension);

            random.setSeed(seed);

            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    randomMatrix.setEntry(i, j, random.nextDouble());
                }
            }
        } while (isSingular(randomMatrix));

        return randomMatrix;
    }

    public static RealMatrix toMatrix(double[] array, int n) {
        RealMatrix matrix = new BlockRealMatrix(n, n);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                matrix.setEntry(i, j, array[i + n * j]);
            }
        }
        return matrix;
    }

    public static void fillArray(RealMatrix matrix, double[] array) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                array[i + matrix.getRowDimension() * j] = matrix.getEntry(i, j);
            }
        }
    }

    public static void fillVector(RealVector vector, double[] array) {
        for (int i = 0; i < vector.getDimension(); i++) {
            array[i] = vector.getEntry(i);
        }
    }

    public static RealVector toVector(double[] array) {
        return new ArrayRealVector(array);
    }

    public static double rescale(double[] array) {
        return Utils.rescale(array, 0);
    }

    public static double rescale(double[] array, double previousLogFactor) {
        double max = Arrays.stream(array).max().orElse(1.0);

        for (int i = 0; i < array.length; i++) {
            array[i] /= max;
        }

        return Math.log(max) + previousLogFactor;
    }
}
