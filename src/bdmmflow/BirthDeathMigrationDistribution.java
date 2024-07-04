package bdmmflow;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.extinctionSystem.ExtinctionProbabilitiesODESystem;
import bdmmflow.flowSystems.*;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.intervals.IntervalUtils;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import beast.base.core.*;
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
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.Arrays;
import java.util.List;

@Citation(value = "Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. (2016). " +
        "A General and Efficient Algorithm for the Likelihood of Diversification and Discrete-Trait Evolutionary Models, \n" +
        "Systematic Biology, Volume 69, Issue 3, May 2020, Pages 545–556."
        , DOI = "10.1093/sysbio/syz055", year = 2020, firstAuthorSurname = "Louca")

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel \" +\n" +
        "        \"with discrete locations and migration events among demes. \" +\n" +
        "        \"This implementation uses the Flow representation of the probability  \" +\n" +
        "        \"ODE for better performance. \" +\n" +
        "        \"It can be used as a drop-in replacement of the BDMM-Prime package."
)
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

    public Input<Boolean> useRandomInitialMatrixInput = new Input<>(
            "useRandomInitialMatrix",
            "Whether to use a random initial matrix as the initial flow state.",
            true
    );

    public Input<Integer> minNumIntervalsInput = new Input<>(
            "minNumIntervals",
            "The number of intervals the time span is broken up into. Increase this when running into numerical stability issues.",
            4
    );

    public Input<Boolean> useInverseFlowInput = new Input<>(
            "useInverseFlow",
            "Whether to use the inverse flow algorithm. It is faster, but can lead to higher numerical instability.",
            true
    );

    public Input<Integer> seedInput = new Input<>(
            "seed",
            "The random seed used in the analysis.",
            3215
    );


    private Parameterization parameterization;

    private boolean useRandomInitialMatrix;

    private double finalSampleOffset;
    private TreeInterface tree;
    private String typeLabel;
    private TraitSet typeTraitSet;

    double[] frequencies;

    boolean conditionOnRoot;
    boolean conditionOnSurvival;

    double absoluteTolerance;
    double relativeTolerance;

    int minNumIntervals;
    boolean useInverseFlow;
    int seed;

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
        this.useRandomInitialMatrix = this.useRandomInitialMatrixInput.get();
        this.seed = this.seedInput.get();
        this.minNumIntervals = this.minNumIntervalsInput.get();
        this.useInverseFlow = this.useInverseFlowInput.get();

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

        // validate minNumIntervals

        if (minNumIntervals < 1) {
            throw new RuntimeException(
                    "Error: minNumIntervals must be at least 1."
            );
        }

        // initialize utils

        this.logScalingFactors = new double[this.tree.getNodeCount()];
        this.initializeIsRhoSampled();
    }

    /**
     * Initializes the isRhoSampled array. The array contains a boolean for every node indicating
     * whether it was rho-sampled or not.
     */
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

    /**
     * Calculates the log tree likelihood.
     *
     * @param dummyTree a dummyTree that is not considered, a BEAST implementation detail.
     * @return the calculated log tree likelihood for the tree specified in the given parameterization.
     */
    @Override
    public double calculateTreeLogLikelihood(TreeInterface dummyTree) {
        // integrate over the extinction probabilities ODE and the flow ODE

        ExtinctionProbabilities extinctionProbabilities = this.calculateExtinctionProbabilities();
        IFlow flow = this.calculateFlow(extinctionProbabilities);

        // recursively traverse the tree to calculate the root likelihood per state

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
            Log.warning("A singular matrix was detected. This is due to numerical instability. Increase minNumIntervals to mitigate this issue.");
            return Double.NEGATIVE_INFINITY;
        }

        // get tree likelihood by a weighted average of the root likelihood per state

        double treeLikelihood = 0.0;
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            treeLikelihood += rootLikelihoodPerState[i] * this.frequencies[i];
        }

        // consider different ways to condition the tree

        double conditionDensity = this.calculateConditionDensityFactor(extinctionProbabilities);
        treeLikelihood /= conditionDensity;

        // turn the likelihood into log likelihood and correct for scaling

        double logTreeLikelihood = Math.log(treeLikelihood) + this.logScalingFactors[root.getNr()];

        // convert from oriented to labeled tree likelihood

        int internalNodeCount = tree.getLeafNodeCount() - ((Tree) tree).getDirectAncestorNodeCount() - 1;
        logTreeLikelihood += Math.log(2) * internalNodeCount - Gamma.logGamma(tree.getLeafNodeCount() + 1);

        return logTreeLikelihood;
    }

    /**
     * Integrates over the extinction probabilities ODE.
     *
     * @return a wrapper class that allows to query the extinction probabilities at any given time.
     */
    private ExtinctionProbabilities calculateExtinctionProbabilities() {
        IntervalODESystem system = new ExtinctionProbabilitiesODESystem(
                this.parameterization,
                this.absoluteTolerance,
                this.relativeTolerance
        );

        // create the initial state

        int endInterval = this.parameterization.getTotalIntervalCount() - 1;

        double[] initialState = new double[this.parameterization.getNTypes()];
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            initialState[i] = 1 - this.parameterization.getRhoValues()[endInterval][i];
        }

        // integrate

        ContinuousOutputModel[] integrationResults = system.integrateBackwards(initialState);

        return new ExtinctionProbabilities(integrationResults);
    }

    /**
     * Precomputes the flow ODE.
     *
     * @param extinctionProbabilities the precomputed extinction probabilities.
     * @return a wrapper class that allows to query the flow at any given time.
     */
    private IFlow calculateFlow(ExtinctionProbabilities extinctionProbabilities) {
        List<Interval> intervals = IntervalUtils.getIntervals(
                this.parameterization,
                this.parameterization.getTotalProcessLength() / this.minNumIntervals
        );

        RealMatrix initialMatrix = this.useRandomInitialMatrix ?
                Utils.getRandomMatrix(this.numTypes, this.seed) :
                MatrixUtils.createRealIdentityMatrix(this.numTypes);
        RealMatrix inverseInitialMatrix = MatrixUtils.inverse(initialMatrix);

        boolean resetInitialStateAtIntervalBoundaries = 1 < this.minNumIntervals;

        IFlowODESystem system;

        if (this.useInverseFlow) {
            system = new InverseFlowODESystem(
                    this.parameterization,
                    extinctionProbabilities,
                    this.absoluteTolerance,
                    this.relativeTolerance
            );
        } else {
            system = new FlowODESystem(
                    this.parameterization,
                    extinctionProbabilities,
                    this.absoluteTolerance,
                    this.relativeTolerance
            );
        }

        return system.calculateFlowIntegral(
                intervals,
                initialMatrix,
                inverseInitialMatrix,
                resetInitialStateAtIntervalBoundaries
        );
    }

    /**
     * Calculates the probability density factor due to the way the tree is conditioned.
     * <p>
     * See Tanja Stadler, How Can We Improve Accuracy of Macroevolutionary Rate Estimates?,
     * Systematic Biology, Volume 62, Issue 2, March 2013, Pages 321–329,
     * https://doi.org/10.1093/sysbio/sys073
     *
     * @param extinctionProbabilities the calculated extinction probabilities integral.
     * @return the factor due to tree conditioning.
     */
    private double calculateConditionDensityFactor(ExtinctionProbabilities extinctionProbabilities) {
        double conditionDensity = 0.0;

        if (this.conditionOnRoot) {
            double[] extinctionAtRoot = extinctionProbabilities.getProbability(0);

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
            double[] extinctionAtRoot = extinctionProbabilities.getProbability(0);

            for (int type = 0; type < parameterization.getNTypes(); type++) {
                conditionDensity += this.frequencies[type] * (1 - extinctionAtRoot[type]);
            }
        } else {
            conditionDensity = 1.0;
        }

        return conditionDensity;
    }

    /**
     * Calculates the per-type likelihood of the subtree of the given node including the edge leading to the
     * node.
     */
    private double[] calculateSubTreeLikelihood(
            Node root,
            double timeEdgeStart,
            double timeEdgeEnd,
            IFlow flow,
            ExtinctionProbabilities extinctionProbabilities
    ) {
        double[] likelihoodEdgeEnd;

        if (root.isLeaf()) {
            likelihoodEdgeEnd = calculateLeafLikelihood(root, timeEdgeEnd, extinctionProbabilities);
        } else if (root.getChild(0).isDirectAncestor() || root.getChild(1).isDirectAncestor()) {
            likelihoodEdgeEnd = calculateDirectAncestorWithChildLikelihood(root, timeEdgeEnd, flow, extinctionProbabilities);
        } else {
            likelihoodEdgeEnd = calculateInternalEdgeLikelihood(root, timeEdgeEnd, flow, extinctionProbabilities);
        }

        return flow.integrateUsingFlow(
                timeEdgeStart,
                timeEdgeEnd,
                likelihoodEdgeEnd
        );
    }

    /**
     * Calculates the likelihood of a single leaf node including the edge leading to it.
     */
    private double[] calculateLeafLikelihood(
            Node root,
            double timeEdgeEnd,
            ExtinctionProbabilities extinctionProbabilities
    ) {

        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);
        double[] extinctionProbabilityEdgeEnd = extinctionProbabilities.getProbability(timeEdgeEnd);

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

    /**
     * Calculates the likelihood of a node that has a direct ancestor as a child; including the edge leading to it.
     */
    private double[] calculateDirectAncestorWithChildLikelihood(
            Node root,
            double timeEdgeEnd,
            IFlow flow,
            ExtinctionProbabilities extinctionProbabilities
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

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        int nodeType = this.getNodeType(directAncestor);

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

    /**
     * Calculates the likelihood of the subtree of the given internal node including the edge leading to it.
     */
    private double[] calculateInternalEdgeLikelihood(
            Node root,
            double timeEdgeEnd,
            IFlow flow,
            ExtinctionProbabilities extinctionProbabilities
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

    /**
     * Returns the type of the given node.
     * @param node the node to return the type of.
     * @return the type of the node.
     */
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

    @Override
    public boolean requiresRecalculation() {
        return true;
    }
}
