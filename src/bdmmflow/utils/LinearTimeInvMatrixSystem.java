package bdmmflow.utils;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;

public class LinearTimeInvMatrixSystem implements FirstOrderDifferentialEquations {

    RealMatrix systemMatrix;
    int matrixDimension;

    public LinearTimeInvMatrixSystem(RealMatrix systemMatrix) {
        this.systemMatrix = systemMatrix;
        this.matrixDimension = systemMatrix.getRowDimension();
    }

    @Override
    public int getDimension() {
        return this.matrixDimension * this.matrixDimension;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        RealMatrix yMatrix = Utils.toMatrix(y, this.matrixDimension);
        RealMatrix yDotMatrix = this.systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }

    public ContinuousOutputModel integrateBackwards(double[] initialState, double startTime, double endTime) {
        ContinuousOutputModel result = new ContinuousOutputModel();

        EulerIntegrator integrator = new EulerIntegrator(
                (endTime - startTime) / 5
        );

        double[] state = initialState.clone();

        integrator.addStepHandler(result);
        integrator.integrate(this, endTime, state, startTime, state);
        integrator.clearStepHandlers();

        return result;
    }
}
