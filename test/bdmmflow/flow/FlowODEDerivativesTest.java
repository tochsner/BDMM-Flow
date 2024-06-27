package bdmmflow.flow;


import bdmmprime.distribution.P0GeSystem;
import bdmmprime.parameterization.*;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;

public class FlowODEDerivativesTest {

    static class NormalFlowODESystem extends FlowODESystem {
        public NormalFlowODESystem(Parameterization parameterization, ContinuousIntervalOutputModel extinctionProbabilities) {
            super(parameterization, extinctionProbabilities, 1e-100, 1e-7);
        }

        @Override
        public int getDimension() {
            return param.getNTypes();
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot) {
            RealMatrix systemMatrix = this.buildSystemMatrix(t);
            RealVector yDotMatrix = systemMatrix.operate(Utils.toVector(y));
            Utils.fillVector(yDotMatrix, yDot);
        }
    }

    void testDerivatives(
            Parameterization parameterization,
            double t,
            double[] initialState
    ) {
        // setup extinction ODE

        IntervalODESystem extinctionSystem = new ExtinctionODESystem(parameterization, 1e-100, 1e-7);

        double[] initialExtinctionState = new double[parameterization.getNTypes()];
        Arrays.fill(initialExtinctionState, 1.0);

        ContinuousIntervalOutputModel extinctionProbabilities = new ContinuousIntervalOutputModel(extinctionSystem.integrateBackwards(
                initialExtinctionState
        ));

        // get flow ODE derivatives

        FlowODESystem flowSystem = new NormalFlowODESystem(parameterization, extinctionProbabilities);
        double[] flowDerivatives = new double[initialState.length];
        flowSystem.computeDerivatives(t, initialState.clone(), flowDerivatives);

        // get old ODE derivatives

        P0GeSystem oldSystem = new P0GeSystem(parameterization, 1e-100, 1e-20);

        double[] startExtinction = extinctionProbabilities.getOutput(t).clone();
        double[] oldInitialState = new double[]{
                startExtinction[0],
                startExtinction[1],
                initialState[0],
                initialState[1]
        };
        double[] oldDerivatives = new double[2*initialState.length];
        oldSystem.computeDerivatives(t, oldInitialState, oldDerivatives);

        for (int i = 0; i < initialState.length; i++) {
            assertEquals(oldDerivatives[i+2], flowDerivatives[i], 1e-7);
        }
    }

    @Test
    public void testFlowODESolver5() {
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("2.5"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(40.0 / 3.0) + " " + Double.toString(1.0 / 3.0)), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(1.0 / 3.0)), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1"),
                        2)
        );

        testDerivatives(
                parameterization,
                1.9044996386612278,
                new double[]{9.127320657394401E-7, 1.655901936494904E-4}
        );
    }

}