package bdmmflow.utils;

import org.apache.commons.math3.linear.*;
import org.jblas.DoubleMatrix;

import java.util.Arrays;
import java.util.Random;

public class Utils {

    /**
     * Returns a random square matrix of the given dimension.
     */
    public static RealMatrix getRandomMatrix(int dimension) {
        return getRandomMatrix(dimension, new Random().nextInt());
    }

    /**
     * Returns a random square matrix of the given dimension.
     */
    public static RealMatrix getRandomMatrix(int dimension, int seed) {
        RealMatrix randomMatrix = new BlockRealMatrix(dimension, dimension);
        Random random = new Random(seed);

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                randomMatrix.setEntry(i, j, random.nextDouble());
            }
        }

        return randomMatrix;
    }

    public static RealMatrix computeRegularMinimizer(RealMatrix A) {
        // Compute SVD
        SingularValueDecomposition svd = new SingularValueDecomposition(A);
        RealMatrix V = svd.getV();
        double[] sigma = svd.getSingularValues();
        int n = sigma.length;

        // Compute product of singular values directly (avoid streams)
        double logProd = 0.0;
        for (double s : sigma) logProd += Math.log(s);
        double geomMean = Math.exp((Math.log(1.0) + logProd) / n); // d=1.0 inline

        // Compute V * Σ⁻¹ efficiently (scale columns of V by 1/σᵢ)
        double[][] vData = V.getData();
        for (int j = 0; j < n; j++) {
            double scale = geomMean / sigma[j];  // combine scalarMultiply(g_d) + Σ⁻¹
            for (int i = 0; i < n; i++) {
                vData[i][j] *= scale;
            }
        }
        RealMatrix Xstar = MatrixUtils.createRealMatrix(vData).multiply(V.transpose());

        return Xstar;
    }

    /**
     * Returns a RealMatrix matrix filled with the values in the array in column-major order.
     */
    public static RealMatrix toMatrix(double[] array, int n) {
        return Utils.toMatrix(array, n, 0);
    }

    /**
     * Returns a RealMatrix matrix filled with the values in the array in column-major order.
     */
    public static RealMatrix toMatrix(double[] array, int n, int offset) {
        RealMatrix matrix = new BlockRealMatrix(n, n);
        double[] column = new double[n];
        for (int j = 0; j < n; j++) {
            System.arraycopy(array, j * n + offset, column, 0, n);
            matrix.setColumn(j, column);
        }
        return matrix;
    }

    /**
     * Converts a packed 1D array of upper triangular elements into a full n x n RealMatrix.
     */
    public static RealMatrix toUpperTriangularMatrix(double[] array, int offset, int n) {
        RealMatrix matrix = new BlockRealMatrix(n, n);
        int index = 0;

        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                // Fill only the upper triangle (j >= i)
                matrix.setEntry(i, j, array[offset + index++]);
            }
            // Note: Lower triangle (i > j) remains 0.0 by default in BlockRealMatrix
        }
        return matrix;
    }

    public static RealMatrix toMatrix(DoubleMatrix source) {
        RealMatrix destination = new BlockRealMatrix(source.rows, source.columns);
        for (int i = 0; i < source.rows; i++) {
            for (int j = 0; j < source.columns; j++) {
                destination.setEntry(i, j, source.get(i, j));
            }
        }
        return destination;
    }

    public static DoubleMatrix toMatrix(RealMatrix source) {
        DoubleMatrix destination = new DoubleMatrix(source.getRowDimension(), source.getColumnDimension());
        for (int i = 0; i < source.getRowDimension(); i++) {
            for (int j = 0; j < source.getColumnDimension(); j++) {
                destination.put(i, j, source.getEntry(i, j));
            }
        }
        return destination;
    }

    /**
     * Fills the given array with the values in the given matrix in column-major order.
     */
    public static void fillArray(RealMatrix matrix, double[] array) {
        Utils.fillArray(matrix, array, 0);
    }

    /**
     * Fills the given array with the values in the given matrix in column-major order.
     */
    public static void fillArray(RealMatrix matrix, double[] array, int offset) {
        int rows = matrix.getRowDimension();
        int cols = matrix.getColumnDimension();
        for (int j = 0; j < cols; j++) {
            System.arraycopy(matrix.getColumn(j), 0, array, j * rows + offset, rows);
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
     *
     * @param array the array to scale.
     * @return the log of the scaling factor applied.
     */
    public static double rescale(double[] array) {
        return Utils.rescale(array, 0);
    }

    /**
     * Scales the given array in-place such that the maximum value is 1. Assumes that
     * the array has already been scaled by the given previousLogFactor.
     *
     * @param array             the array to scale.
     * @param previousLogFactor the log of the previous scaling factor.
     * @return the log of the overall scaling factor applied.
     */
    public static double rescale(double[] array, double previousLogFactor) {
        double max = Arrays.stream(array).max().orElse(1.0);

        if (max == 0.0) return previousLogFactor;

        for (int i = 0; i < array.length; i++) {
            array[i] /= max;
        }

        return Math.log(max) + previousLogFactor;
    }
}
