//package bdmmflow.flow;
//
//
//import bdmmflow.flow.extinctionSystem.ExtinctionProbabilities;
//import bdmmflow.flow.extinctionSystem.ExtinctionProbabilitiesODESystem;
//import bdmmflow.flow.flowSystems.InverseFlowODESystem;
//import bdmmflow.flow.intervals.IntervalODESystem;
//import bdmmflow.flow.utils.Utils;
//import bdmmprime.distribution.P0GeSystem;
//import bdmmprime.parameterization.*;
//import beast.base.inference.parameter.RealParameter;
//import org.apache.commons.math3.linear.RealMatrix;
//import org.apache.commons.math3.linear.RealVector;
//import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
//import org.apache.commons.math3.ode.FirstOrderIntegrator;
//import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
//import org.junit.Test;
//
//import java.util.Arrays;
//
//import static junit.framework.Assert.assertEquals;
//import static org.junit.jupiter.api.Assertions.assertArrayEquals;
//
//public class FlowODESystemIsolatedTest {
//
//    static class NormalFlowODESystem extends InverseFlowODESystem {
//        public NormalFlowODESystem(Parameterization parameterization, ExtinctionProbabilities extinctionProbabilities) {
//            super(parameterization, extinctionProbabilities, 1e-100, 1e-20);
//        }
//
//        @Override
//        public int getDimension() {
//            return param.getNTypes();
//        }
//
//        @Override
//        public void computeDerivatives(double t, double[] y, double[] yDot) {
//            RealMatrix systemMatrix = this.buildSystemMatrix(t);
//            RealVector yDotMatrix = systemMatrix.operate(Utils.toVector(y));
//            Utils.fillVector(yDotMatrix, yDot);
//        }
//    }
//
//    public void testFlowODESolver(
//            Parameterization parameterization,
//            double[] initialState,
//            double startTime,
//            double endTime
//    ) {
//        // set up extinction ODE
//
//        IntervalODESystem extinctionSystem = new ExtinctionProbabilitiesODESystem(parameterization, 1e-100, 1e-20);
//
//        double[] initialExtinctionState = new double[parameterization.getNTypes()];
//        Arrays.fill(initialExtinctionState, 1.0);
//
//        ExtinctionProbabilities extinctionProbabilities = new ExtinctionProbabilities(extinctionSystem.integrateBackwards(
//                initialExtinctionState
//        ));
//
//        int intervalEndTime = parameterization.getIntervalIndex(endTime);
//        double[] endExtinction = extinctionProbabilities.getProbability(endTime).clone();
//
//        // set up integrator
//
//        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
//                parameterization.getTotalProcessLength() * 1e-100,
//                parameterization.getTotalProcessLength() / 10, 1e-100, 1e-10
//        );
//
//        // integrate using flow
//
//        FirstOrderDifferentialEquations flowSystem = new NormalFlowODESystem(
//                parameterization, extinctionProbabilities
//        );
//        double[] flowIntegral = new double[2];
//        integrator.integrate(flowSystem, endTime, initialState.clone(), startTime, flowIntegral);
//
//        // integrate using normal method
//
//        FirstOrderDifferentialEquations normalSystem = new P0GeSystem(
//                parameterization, 1e-100, 1e-10
//        );
//        double[] normalIntegral = new double[4];
//        double[] normalInitialState = new double[] {
//                endExtinction[0], endExtinction[1], initialState[0], initialState[1]
//        };
//        integrator.integrate(normalSystem, endTime, normalInitialState.clone(), startTime, normalIntegral);
//        normalIntegral = Arrays.copyOfRange(normalIntegral, 2, 4);
//
//        // assert equality
//
//        System.out.println(Arrays.toString(normalIntegral));
//        System.out.println(Arrays.toString(flowIntegral));
//        assertArrayEquals(normalIntegral, flowIntegral, 1e-10);
//    }
//
//    @Test
//    public void testFlowODE() {
//        Parameterization parameterization = new EpiParameterization();
//        parameterization.initByName(
//                "processLength", new RealParameter("2.5"),
//                "typeSet", new TypeSet(2),
//                "R0", new SkylineVectorParameter(
//                        null,
//                        new RealParameter(Double.toString(40.0 / 3.0) + " " + Double.toString(0.0 / 3.0)), 2),
//                "becomeUninfectiousRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.5"), 2),
//                "samplingProportion", new SkylineVectorParameter(
//                        null,
//                        new RealParameter(Double.toString(1.0 / 3.0)), 2),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.0"), 2),
//                "migrationRate", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.1"),
//                        2)
//        );
//
//        double[] initialState = new double[]{0.5, 0};
//
//        testFlowODESolver(parameterization, initialState, 1.0, 2.5);
//    }
//
//    @Test
//    public void testFlowODE2() {
//        Parameterization parameterization = new EpiParameterization();
//        parameterization.initByName(
//                "typeSet", new TypeSet(2),
//                "processLength", new RealParameter("6.0"),
//                "R0", new SkylineVectorParameter(
//                        null,
//                        new RealParameter((4.0/3.0) + " 1.1"),
//                        2),
//                "becomeUninfectiousRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.5 1.4"),
//                        2),
//                "R0AmongDemes", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.0"),
//                        2),
//                "migrationRate", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.2 0.3"),
//                        2),
//                "samplingProportion", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.33"),
//                        2),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.3 0.4"))
//        );
//
//        double[] initialState = new double[]{0E0, 11.35022270351E-3};
//
//        testFlowODESolver(parameterization, initialState, 1.0, 4.5);
//    }
//
//}