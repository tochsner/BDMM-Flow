package bdmmflow.intervals;

import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalUtils;
import bdmmprime.parameterization.*;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;
import org.junit.jupiter.api.Assertions;

import java.util.Arrays;
import java.util.List;

public class GetIntervalsTest {
    @Test
    public void testOnlyParameterizationIntervals() {
        Parameterization parameterization = new CanonicalParameterization();

        RealParameter changeTimes = new RealParameter("1.0 3.0 7.0 9.0");

        parameterization.initByName(
                "processLength", new RealParameter("10.0"),
                "birthRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("2.0 2.0 2.0 2.0 2.0")
                ),
                "deathRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                ),
                "samplingRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                ),
                "removalProb", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                )
        );

        List<Interval> actualIntervals = IntervalUtils.getIntervals(parameterization, 10);
        List<Interval> expectedIntervals = List.of(
                new Interval(0, List.of(0, 1, 2, 3), 0.0, 10.0)
        );

        Assertions.assertIterableEquals(expectedIntervals, actualIntervals);
    }

    @Test
    public void testOnlySmallEnoughIntervals() {
        Parameterization parameterization = new CanonicalParameterization();

        RealParameter changeTimes = null;

        parameterization.initByName(
                "processLength", new RealParameter("10.0"),
                "birthRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("2.0")
                ),
                "deathRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5")
                ),
                "samplingRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5")
                ),
                "removalProb", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5")
                )
        );

        List<Interval> actualIntervals = IntervalUtils.getIntervals(parameterization, 2);
        List<Interval> expectedIntervals = Arrays.asList(
                new Interval(0, List.of(), 0.0, 2.0),
                new Interval(1, List.of(), 2.0, 4.0),
                new Interval(2, List.of(), 4.0, 6.0),
                new Interval(3, List.of(), 6.0, 8.0),
                new Interval(4, List.of(), 8.0, 10.0)
        );

        Assertions.assertIterableEquals(expectedIntervals, actualIntervals);
    }

    @Test
    public void testCombination() {
        Parameterization parameterization = new CanonicalParameterization();

        RealParameter changeTimes = new RealParameter("1.0 3.0 7.0 9.0");

        parameterization.initByName(
                "processLength", new RealParameter("10.0"),
                "birthRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("2.0 2.0 2.0 2.0 2.0")
                ),
                "deathRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                ),
                "samplingRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                ),
                "removalProb", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                )
        );

        List<Interval> actualIntervals = IntervalUtils.getIntervals(parameterization, 2);
        List<Interval> expectedIntervals = Arrays.asList(
                new Interval(0, List.of(0), 0.0, 2.0),
                new Interval(1, List.of(1), 2.0, 4.0),
                new Interval(2, List.of(), 4.0, 6.0),
                new Interval(3, List.of(2), 6.0, 8.0),
                new Interval(4, List.of(3), 8.0, 10.0)
        );

        Assertions.assertIterableEquals(expectedIntervals, actualIntervals);
    }

    @Test
    public void testCombinationWithAlignedIntervals() {
        Parameterization parameterization = new CanonicalParameterization();

        RealParameter changeTimes = new RealParameter("2.0 3.0 7.0 9.0");

        parameterization.initByName(
                "processLength", new RealParameter("10.0"),
                "birthRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("2.0 2.0 2.0 2.0 2.0")
                ),
                "deathRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                ),
                "samplingRate", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                ),
                "removalProb", new SkylineVectorParameter(
                        changeTimes,
                        new RealParameter("0.5 0.5 0.5 0.5 0.5")
                )
        );

        List<Interval> actualIntervals = IntervalUtils.getIntervals(parameterization, 2);
        List<Interval> expectedIntervals = Arrays.asList(
                new Interval(0, List.of(), 0.0, 2.0),
                new Interval(1, List.of(1), 2.0, 4.0),
                new Interval(2, List.of(), 4.0, 6.0),
                new Interval(3, List.of(2), 6.0, 8.0),
                new Interval(4, List.of(3), 8.0, 10.0)
        );

        Assertions.assertIterableEquals(expectedIntervals, actualIntervals);
    }
}