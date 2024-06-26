package bdmmprime.flow.initialMatrices;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import java.util.Random;


public class InitialMatrices {

    public enum MatrixType {
        IDENTITY,
        RANDOM,
        GAUSSIAN,
        EXP
    }

    static Random random = new Random();

    static boolean isSingular(RealMatrix matrix) {
        if (matrix.getRowDimension() != matrix.getColumnDimension()) return false;

        SingularValueDecomposition svd = new SingularValueDecomposition(matrix);
        return svd.getRank() != matrix.getColumnDimension();
    }

    public static RealMatrix getIdentity(int dimension) {
        return MatrixUtils.createRealIdentityMatrix(dimension);
    }

    public static RealMatrix getRandom(int dimension) {
        RealMatrix randomMatrix;

        do {
            randomMatrix = new BlockRealMatrix(dimension, dimension);

            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    randomMatrix.setEntry(i, j, random.nextDouble());
                }
            }
        } while (isSingular(randomMatrix));

        return randomMatrix;
    }

    public static RealMatrix getGaussian(int dimension) {
        RealMatrix randomMatrix;

        do {
            randomMatrix = new BlockRealMatrix(dimension, dimension);

            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    randomMatrix.setEntry(i, j, random.nextGaussian());
                }
            }
        } while (isSingular(randomMatrix));

        return randomMatrix;
    }

    public static RealMatrix getExp(int dimension) {
        RealMatrix randomMatrix;

        do {
            randomMatrix = new BlockRealMatrix(dimension, dimension);

            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    randomMatrix.setEntry(i, j, random.nextExponential());
                }
            }
        } while (isSingular(randomMatrix));

        return randomMatrix;
    }

    public static RealMatrix getMatrix(int dimension, MatrixType type) {
        return switch (type) {
            case GAUSSIAN -> InitialMatrices.getGaussian(dimension);
            case EXP -> InitialMatrices.getExp(dimension);
            case RANDOM -> InitialMatrices.getRandom(dimension);
            default -> InitialMatrices.getIdentity(dimension);
        };
    }
}
