package bdmmflow.utils;

import org.apache.commons.math3.linear.*;

import java.util.Arrays;
import java.util.Random;

public class Utils {

    private static final Random random = new Random();

    /**
     * Returns whether the matrix is singular.
     * @param matrix the matrix to check.
     * @return if the matrix is singular.
     */
    public static boolean isSingular(RealMatrix matrix) {
        if (matrix.getRowDimension() != matrix.getColumnDimension()) return false;

        SingularValueDecomposition svd = new SingularValueDecomposition(matrix);
        return svd.getRank() != matrix.getColumnDimension();
    }

    /**
     * Returns a random square matrix of the given dimension.
     */
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

    /**
     * Returns a RealMatrix matrix filled with the values in the array in column-major order.
     */
    public static RealMatrix toMatrix(double[] array, int n) {
        RealMatrix matrix = new BlockRealMatrix(n, n);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                matrix.setEntry(i, j, array[i + n * j]);
            }
        }
        return matrix;
    }

    /**
     * Fills the given array with the values in the given matrix in column-major order.
     */
    public static void fillArray(RealMatrix matrix, double[] array) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                array[i + matrix.getRowDimension() * j] = matrix.getEntry(i, j);
            }
        }
    }

    /**
     * Converts the given array to a RealVector vector.
     */
    public static RealVector toVector(double[] array) {
        return new ArrayRealVector(array);
    }

    /**
     * Scales the given array in-place such that the maximum value is 1.
     * @param array the array to scale.
     * @return the log of the scaling factor applied.
     */
    public static double rescale(double[] array) {
        return Utils.rescale(array, 0);
    }

    /**
     * Scales the given array in-place such that the maximum value is 1. Assumes that
     * the array has already been scaled by the given previousLogFactor.
     * @param array the array to scale.
     * @param previousLogFactor the log of the previous scaling factor.
     * @return the log of the overall scaling factor applied.
     */
    public static double rescale(double[] array, double previousLogFactor) {
        double max = Arrays.stream(array).max().orElse(1.0);

        for (int i = 0; i < array.length; i++) {
            // array[i] /= max;
        }

        return 0; // Math.log(max) + previousLogFactor;
    }
}
