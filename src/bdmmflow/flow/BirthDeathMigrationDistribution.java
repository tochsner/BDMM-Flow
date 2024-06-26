package bdmmflow.flow;

import bdmmprime.flow.initialMatrices.InitialMatrices;
import beast.base.core.Log;
import bdmmprime.parameterization.Parameterization;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput = new Input<>(
            "parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED
    );

    public Input<Function> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0")
    );

    public Input<RealParameter> frequenciesInput = new Input<>(
            "frequencies",
            "The equilibrium frequencies for each type",
            new RealParameter("1.0")
    );

    public Input<String> typeLabelInput = new Input<>(
            "typeLabel",
            "Attribute key used to specify sample trait values in tree."
    );

    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet",
            "Trait set specifying sample trait values.");

    public Input<Boolean> conditionOnSurvivalInput = new Input<>("conditionOnSurvival",
            "Condition on at least one surviving lineage. (Default true.)",
            true);

    public Input<Boolean> conditionOnRootInput = new Input<>("conditionOnRoot",
            "Condition on root age, not time of origin.", false);

    public Input<Double> relativeToleranceInput = new Input<>(
            "relTolerance",
            "Relative tolerance for numerical integration.",
            1e-7
    );

    public Input<Double> absoluteToleranceInput = new Input<>(
            "absTolerance",
            "Absolute tolerance for numerical integration.",
            1e-100
    );

    /**
     * DEBUG STIFFNESS
     */

    public Input<Boolean> useIntervalsInput = new Input<>(
            "useIntervals",
            "",
            true
    );
    boolean useIntervals;

    public Input<Integer> minNumIntervalsInput = new Input<>(
            "minNumIntervals",
            "",
            10
    );
    int minNumIntervals;

    public Input<String> integratorInput = new Input<>(
            "integrator",
            "",
            "DormandPrince54Integrator"
    );
    String integrator;

    public Input<String> initialMatrixTypeInput = new Input<>(
            "initialMatrixType",
            "",
            "RANDOM"
    );
    InitialMatrices.MatrixType initialMatrixType;


    private Parameterization parameterization;

    private double finalSampleOffset;
    private TreeInterface tree;
    private String typeLabel;
    private TraitSet typeTraitSet;

    double[] frequencies;

    boolean conditionOnRoot;
    boolean conditionOnSurvival;

    double absoluteTolerance;
    double relativeTolerance;

    int numTypes;

    double[] logScalingFactors;
    boolean[] isRhoSampled;

    @Override
    public void initAndValidate() {
        // unpack input values

        this.parameterization = this.parameterizationInput.get();
        this.finalSampleOffset = this.finalSampleOffsetInput.get().getArrayValue();
        this.tree = this.treeInput.get();
        this.typeLabel = this.typeLabelInput.get();
        this.typeTraitSet = this.typeTraitSetInput.get();
        this.frequencies = this.frequenciesInput.get().getDoubleValues();
        this.conditionOnRoot = this.conditionOnRootInput.get();
        this.conditionOnSurvival = this.conditionOnSurvivalInput.get();
        this.absoluteTolerance = this.absoluteToleranceInput.get();
        this.relativeTolerance = this.relativeToleranceInput.get();
        this.numTypes = this.parameterization.getNTypes();

        this.minNumIntervals = this.minNumIntervalsInput.get();
        this.useIntervals = this.useIntervalsInput.get();
        this.integrator = this.integratorInput.get();
        this.initialMatrixType = InitialMatrices.MatrixType.valueOf(this.initialMatrixTypeInput.get());

        // validate typeLabelInput

        if (this.numTypes != 1 && this.typeLabel == null && this.typeTraitSet == null) {
            throw new RuntimeException(
                    "Error: For models with >1 type, typeLabel or typeTraitSet must be specified."
            );
        }

        // validate frequenciesInput

        if (this.frequencies.length != numTypes) {
            throw new RuntimeException(
                    "Error: dimension of equilibrium frequencies parameter must match number of types."
            );
        }

        double freqSum = 0.0;
        for (double f : this.frequencies) {
            freqSum += f;
        }
        if (!bdmmprime.util.Utils.equalWithPrecision(freqSum, 1.0)) {
            throw new RuntimeException(
                    "Error: equilibrium frequencies must add up to 1 but currently add to %f.".formatted(freqSum)
            );
        }

        // initialize utils needed

        this.logScalingFactors = new double[this.tree.getNodeCount()];
        this.initializeIsRhoSampled();
    }

    private void initializeIsRhoSampled() {
        this.isRhoSampled = new boolean[this.tree.getLeafNodeCount()];

        for (int i = 0; i < this.tree.getLeafNodeCount(); i++) {
            this.isRhoSampled[i] = false;

            double nodeTime = this.parameterization.getNodeTime(this.tree.getNode(i), this.finalSampleOffset);

            for (double rhoSamplingTime : this.parameterization.getRhoSamplingTimes()) {
                if (bdmmprime.util.Utils.equalWithPrecision(nodeTime, rhoSamplingTime)) {
                    this.isRhoSampled[i] = true;
                    break;
                }
            }
        }
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface dummyTree) {
        ContinuousIntervalOutputModel extinctionProbabilities = this.calculateExtinctionProbabilities();
        Flow flow = this.calculateFlow(extinctionProbabilities);

        Node root = this.tree.getRoot();
        double[] rootLikelihoodPerState;
        try {
            rootLikelihoodPerState = this.calculateSubTreeLikelihood(
                    root,
                    0,
                    this.parameterization.getNodeTime(root, this.finalSampleOffset),
                    flow,
                    extinctionProbabilities
            );
        } catch (SingularMatrixException exception) {
            Log.warning("Singular matrix was found during tree likelihood calculation.");
            return Float.NEGATIVE_INFINITY;
        }

        // get tree likelihood by a weighted average of the root likelihood per state

        double treeLikelihood = 0.0;
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            treeLikelihood += rootLikelihoodPerState[i] * this.frequencies[i];
        }

        // consider different ways to condition the tree

        double conditionDensity = this.calculateConditionDensity(extinctionProbabilities);
        treeLikelihood /= conditionDensity;

        double logTreeLikelihood = Math.log(treeLikelihood) + this.logScalingFactors[root.getNr()];

        // convert from oriented to labeled tree likelihood

        int internalNodeCount = tree.getLeafNodeCount() - ((Tree) tree).getDirectAncestorNodeCount() - 1;
        logTreeLikelihood += Math.log(2) * internalNodeCount - Gamma.logGamma(tree.getLeafNodeCount() + 1);

        return logTreeLikelihood;
    }

    private double calculateConditionDensity(ContinuousIntervalOutputModel extinctionProbabilities) {
        // see Tanja Stadler, How Can We Improve Accuracy of Macroevolutionary Rate Estimates?,
        // Systematic Biology, Volume 62, Issue 2, March 2013, Pages 321–329,
        // https://doi.org/10.1093/sysbio/sys073
        double conditionDensity = 0.0;

        if (this.conditionOnRoot) {
            double[] extinctionAtRoot = extinctionProbabilities.getOutput(0);

            int startInterval = this.parameterization.getIntervalIndex(0);

            for (int type1 = 0; type1 < parameterization.getNTypes(); type1++) {
                for (int type2 = 0; type2 < parameterization.getNTypes(); type2++) {
                    double rate = type1 == type2
                            ? parameterization.getBirthRates()[startInterval][type1]
                            : parameterization.getCrossBirthRates()[startInterval][type1][type2];

                    conditionDensity += rate * this.frequencies[type1]
                            * (1 - extinctionAtRoot[type1])
                            * (1 - extinctionAtRoot[type2]);
                }
            }
        } else if (this.conditionOnSurvival) {
            double[] extinctionAtRoot = extinctionProbabilities.getOutput(0);

            for (int type = 0; type < parameterization.getNTypes(); type++) {
                conditionDensity += this.frequencies[type] * (1 - extinctionAtRoot[type]);
            }
        } else {
            conditionDensity = 1.0;
        }

        return conditionDensity;
    }

    private ContinuousIntervalOutputModel calculateExtinctionProbabilities() {
        IntervalODESystem system = new ExtinctionODESystem(this.parameterization, this.absoluteTolerance, this.relativeTolerance);

        int endInterval = this.parameterization.getIntervalIndex(this.parameterization.getTotalProcessLength());

        double[] initialState = new double[this.parameterization.getNTypes()];
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            initialState[i] = 1 - this.parameterization.getRhoValues()[endInterval][i];
        }

        return new ContinuousIntervalOutputModel(system.integrateBackwardsOverIntegrals(initialState));
    }

    private Flow calculateFlow(ContinuousIntervalOutputModel extinctionProbabilities) {
        FlowODESystem system = new FlowODESystem(this.parameterization, extinctionProbabilities, this.absoluteTolerance, this.relativeTolerance, this.integrator);

        RealMatrix initialMatrix = InitialMatrices.getMatrix(this.parameterization.getNTypes(), this.initialMatrixType);
        RealMatrix inverseInitialMatrix = MatrixUtils.inverse(initialMatrix);

        double[] initialState = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
        double[] inverseInitialState = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];

        Utils.fillArray(initialMatrix, initialState);
        Utils.fillArray(inverseInitialMatrix, inverseInitialState);

        return new Flow(
                system.integrateOverIntegrals(
                        initialState, this.useIntervals ? this.parameterization.getTotalProcessLength() / this.minNumIntervals : this.parameterization.getTotalProcessLength(), this.useIntervals
                ), this.numTypes, Utils.toMatrix(inverseInitialState, this.parameterization.getNTypes()), this.useIntervals
        );
    }

    private double[] calculateSubTreeLikelihood(
            Node root,
            double timeEdgeStart,
            double timeEdgeEnd,
            Flow flow,
            ContinuousIntervalOutputModel extinctionProbabilities
    ) {
        double[] likelihoodEdgeEnd;

        if (root.isLeaf()) {
            likelihoodEdgeEnd = calculateLeafLikelihood(root, timeEdgeEnd, extinctionProbabilities);
        } else if (root.getChild(0).isDirectAncestor() || root.getChild(1).isDirectAncestor()) {
            likelihoodEdgeEnd = calculateDirectAncestorWithChildLikelihood(root, timeEdgeEnd, flow, extinctionProbabilities);
        } else {
            likelihoodEdgeEnd = calculateInternalEdgeLikelihood(root, timeEdgeEnd, flow, extinctionProbabilities);
        }

        int intervalEdgeStart = this.parameterization.getIntervalIndex(timeEdgeStart);
        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);

        if (intervalEdgeStart < this.parameterization.getTotalIntervalCount() - 1 && timeEdgeStart == this.parameterization.getIntervalEndTimes()[intervalEdgeStart]) {
            intervalEdgeStart++;
        }

        return FlowODESystem.integrateUsingFlow(
                timeEdgeStart,
                intervalEdgeStart,
                timeEdgeEnd,
                intervalEdgeEnd,
                likelihoodEdgeEnd,
                flow
        );
    }

    private double[] calculateLeafLikelihood(
            Node root,
            double timeEdgeEnd,
            ContinuousIntervalOutputModel extinctionProbabilities
    ) {

        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);
        double[] extinctionProbabilityEdgeEnd = extinctionProbabilities.getOutput(timeEdgeEnd);

        int nodeType = this.getNodeType(root);

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        if (isRhoSampled[root.getNr()]) {
            likelihoodEdgeEnd[nodeType] = this.parameterization.getRhoValues()[intervalEdgeEnd][nodeType];
        } else {
            likelihoodEdgeEnd[nodeType] = this.parameterization.getSamplingRates()[intervalEdgeEnd][nodeType] *
                    (this.parameterization.getRemovalProbs()[intervalEdgeEnd][nodeType]
                            + (1 - this.parameterization.getRemovalProbs()[intervalEdgeEnd][nodeType])
                            * extinctionProbabilityEdgeEnd[nodeType]);
        }

        this.logScalingFactors[root.getNr()] = Utils.rescale(likelihoodEdgeEnd);

        return likelihoodEdgeEnd;
    }

    private double[] calculateDirectAncestorWithChildLikelihood(
            Node root,
            double timeEdgeEnd,
            Flow flow,
            ContinuousIntervalOutputModel extinctionProbabilities
    ) {
        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);

        Node directAncestor = root.getChild(0).isDirectAncestor() ?
                root.getChild(0) : root.getChild(1);
        Node child = root.getChild(0).isDirectAncestor() ?
                root.getChild(1) : root.getChild(0);

        double[] likelihoodChild = this.calculateSubTreeLikelihood(
                child,
                timeEdgeEnd,
                this.parameterization.getNodeTime(child, this.finalSampleOffset),
                flow,
                extinctionProbabilities
        );

        int nodeType = this.getNodeType(directAncestor);

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        if (isRhoSampled[directAncestor.getNr()]) {
            likelihoodEdgeEnd[nodeType] = this.parameterization.getRhoValues()[intervalEdgeEnd][nodeType];
        } else {
            likelihoodEdgeEnd[nodeType] = this.parameterization.getSamplingRates()[intervalEdgeEnd][nodeType];
        }

        likelihoodEdgeEnd[nodeType] *= (1 - this.parameterization.getRemovalProbs()[intervalEdgeEnd][nodeType])
                * likelihoodChild[nodeType];

        this.logScalingFactors[root.getNr()] = Utils.rescale(likelihoodEdgeEnd, this.logScalingFactors[child.getNr()]);

        return likelihoodEdgeEnd;
    }

    private double[] calculateInternalEdgeLikelihood(
            Node root,
            double timeEdgeEnd,
            Flow flow,
            ContinuousIntervalOutputModel extinctionProbabilities
    ) {
        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);

        Node child1 = root.getChild(0);
        double[] likelihoodChild1 = this.calculateSubTreeLikelihood(
                child1,
                timeEdgeEnd,
                this.parameterization.getNodeTime(child1, this.finalSampleOffset),
                flow,
                extinctionProbabilities
        );

        Node child2 = root.getChild(1);
        double[] likelihoodChild2 = this.calculateSubTreeLikelihood(
                child2,
                timeEdgeEnd,
                this.parameterization.getNodeTime(child2, this.finalSampleOffset),
                flow,
                extinctionProbabilities
        );

        // combine the child likelihoods to get the likelihood at the edge end

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            likelihoodEdgeEnd[i] += this.parameterization.getBirthRates()[intervalEdgeEnd][i] * (
                    likelihoodChild1[i] * likelihoodChild2[i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
                if (i == j) {
                    continue;
                }

                likelihoodEdgeEnd[i] += 0.5 * this.parameterization.getCrossBirthRates()[intervalEdgeEnd][i][j] * (
                        likelihoodChild1[i] * likelihoodChild2[j] + likelihoodChild1[j] * likelihoodChild2[i]
                );
            }
        }

        this.logScalingFactors[root.getNr()] = Utils.rescale(
                likelihoodEdgeEnd,
                this.logScalingFactors[child1.getNr()] + this.logScalingFactors[child2.getNr()]
        );

        return likelihoodEdgeEnd;
    }

    private int getNodeType(Node node) {
        if (parameterization.getNTypes() == 1) {
            return 0;
        }

        String nodeTypeName;

        if (this.typeTraitSet != null)
            nodeTypeName = this.typeTraitSet.getStringValue(node.getID());
        else {
            Object metaData = node.getMetaData(this.typeLabel);
            if (metaData instanceof Double) {
                nodeTypeName = String.valueOf(Math.round((double) metaData));
            } else {
                nodeTypeName = metaData.toString();
            }
        }

        return parameterization.getTypeSet().getTypeIndex(nodeTypeName);
    }
}
