from unittest import TestCase

from civaclib.common import TESTDATA, titerSteps, VIRUSES, SERA, LIMIT
from civaclib.parseTiters import (
    getPRNTDiscrete,
    parseTiterExcel,
    averagePlaqueCounts,
    calculateFractionInfectivity,
    convertRawCountsToNeutcurveDf,
    convertRawCountsToDiscreteDf,
    getPRNTContinuous,
)


class TestGetPRNTDiscrete(TestCase):
    """
    Tests for the getPRNTDiscrete function.
    """

    def testFullCurve(self):
        """
        If a full titer curve is given, the correct PRNT50 titer must be read
        out.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 100,
            "1:320": 95.55555556,
            "1:640": 66.22222222,
            "1:1280": 17.33333333,
            "1:2560": 0,
            "1:5120": 0,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("640", prnt50)

    def testGreaterThan(self):
        """
        If the reduction in neutralisation at the highest dilution is still
        above the given threshold, the PRNT50 titer must be expressed as '>'
        the highest dilution.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 100,
            "1:320": 100,
            "1:640": 100,
            "1:1280": 100,
            "1:2560": 100,
            "1:5120": 88.73239437,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual(">5120", prnt50)

    def testLessThanAtThreshold(self):
        """
        If even at the lowest dilution, the reduction in neutralisation
        doesn't rise above the threshold, the titer should be expressed as
        '<' the lowest dilution.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 31.55555556,
            "1:320": 0,
            "1:640": 0,
            "1:1280": 0,
            "1:2560": 0,
            "1:5120": 0,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<160")
        self.assertEqual("<160", prnt50)

    def testLessThanHigherThreshold(self):
        """
        Correct titer must be calculated if it's a less than titer that's
        higher than the specified threshold titer.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 31.55555556,
            "1:320": 0,
            "1:640": 0,
            "1:1280": 0,
            "1:2560": 0,
            "1:5120": 0,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("<160", prnt50)

    def testRoundingUp(self):
        """
        Rounding up to nearest dilution must work.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 100,
            "1:320": 80,
            "1:640": 60,
            "1:1280": 45,
            "1:2560": 20,
            "1:5120": 10,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("1280", prnt50)

    def testRoundingDown(self):
        """
        Rounding down to nearest dilution must work.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 100,
            "1:320": 80,
            "1:640": 52,
            "1:1280": 45,
            "1:2560": 20,
            "1:5120": 10,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("640", prnt50)

    def testHighestDilutionBelowLimit(self):
        """
        If the first dilution is already below the limit, set accordingly.
        """
        # Test this with and without the limit of 5. > than titer
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 46,
            "1:320": 30,
            "1:640": 20,
            "1:1280": 15,
            "1:2560": 12,
            "1:5120": 10,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("160", prnt50)

        data2 = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 44,
            "1:320": 30,
            "1:640": 20,
            "1:1280": 15,
            "1:2560": 12,
            "1:5120": 10,
        }

        prnt502 = getPRNTDiscrete(data2, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("<160", prnt502)

    def testLessThanTiter(self):
        """
        Neutralization threshold is never reached. Treat accordingly.
        """
        data = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 100,
            "1:320": 100,
            "1:640": 100,
            "1:1280": 99,
            "1:2560": 88,
            "1:5120": 77,
        }

        prnt50 = getPRNTDiscrete(data, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual(">5120", prnt50)

        data2 = {
            "1:20": "nd",
            "1:40": "nd",
            "1:80": "nd",
            "1:160": 100,
            "1:320": 100,
            "1:640": 100,
            "1:1280": 99,
            "1:2560": 88,
            "1:5120": 54,
        }

        prnt502 = getPRNTDiscrete(data2, limit=LIMIT, steps=titerSteps, nd="<20")
        self.assertEqual("5120", prnt502)


class TestParseTiterExcel(TestCase):
    """
    Test for the parseTiterExcel function.
    """

    def testCorrectData(self):
        """
        The excel file must be parsed correctly.
        """
        data = parseTiterExcel(TESTDATA, 10, 3, VIRUSES, addAverage=True)

        # Check that the dataframe has the correct columns
        columns = [
            "sample ID",
            "Infektion Hamster",
            "Antigen",
            "Viruskontrolle",
            "Replicate",
            "1:20",
            "1:40",
            "1:80",
            "1:160",
            "1:320",
            "1:640",
            "1:1280",
            "1:2560",
            "1:5120",
        ]
        self.assertEqual(columns, list(data.columns))
        self.assertEqual(len(columns), data.shape[1])
        # Check that the dataframe has the correct length (for each serum
        # and virus pair there are three lines (repeat 1, repeat 2, average)).
        self.assertEqual(len(set(data["sample ID"])) * 3 * len(VIRUSES), data.shape[0])
        # Check that the correct viruses are present in the dataframe
        self.assertEqual(set(VIRUSES), set(data["Antigen"]))

        # Check that the correct sera are present in the dataframe
        self.assertEqual(
            set(SERA),
            set(data["sample ID"]),
        )

        # Check some specific cases
        self.assertEqual(
            [
                "1.1",
                "BA.2 Isolate 26729_2",
                "SARS-CoV-2_WT (984)",
                56.25,
                "PFU Ansatz 1",
                "nd",
                "nd",
                "nd",
                "36",
                ">50",
                ">50",
                ">50",
                ">50",
                ">50",
            ],
            list(data.loc[0]),
        )

        self.assertEqual(
            [
                "6.3",
                "BA.2 Isolate 26729_12",
                "SARS-CoV-2_BA.5",
                111.83,
                "PFU Ansatz 1",
                "nd",
                "nd",
                "nd",
                "16",
                "42",
                "69",
                "103",
                "111",
                "111",
            ],
            list(data.loc[554]),
        )

        self.assertEqual(
            [
                "6.2",
                "BA.2 Isolate 26729_12",
                "SARS-CoV-2_Mu",
                45.83,
                "PFU Ansatz 1",
                "nd",
                "nd",
                "0",
                "2",
                "9",
                "19",
                "35",
                "41",
                "46",
            ],
            list(data.loc[610]),
        )


class TestAveragePlaqueCounts(TestCase):
    """
    Tests for the averagePlaqueCounts function.
    """

    def testAllnd(self):
        """
        Correct result if all values are 'nd'.
        """
        self.assertEqual("nd", averagePlaqueCounts(["nd", "nd"], 53))

    def testOnendOneNumber(self):
        """
        Correct result if one value is nd, one's a number.
        """
        self.assertEqual(40, averagePlaqueCounts(["nd", "40"], 53))

    def testgreaterThan(self):
        """
        Correct result if plaque count is a greater than number.
        """
        self.assertEqual(">50", averagePlaqueCounts([">50", ">50"], 53))

    def testGreaterThanNd(self):
        """
        Correct result if one value is nd, one's a >.
        """
        self.assertEqual(">50", averagePlaqueCounts([">50", "nd"], 53))

    def testAllNumbers(self):
        """
        Correct result if all values are actual numbers.
        """
        self.assertEqual(25, averagePlaqueCounts(["20", "30"], 53))

    def testAllNumbersInt(self):
        """
        Correct result if all values are actual numbers and given as ints.
        """
        self.assertEqual(25, averagePlaqueCounts([20, 30], 53))

    def testgreaterThanNumber(self):
        """
        Correct result if one value is a number, one's a >.
        """
        self.assertEqual(41.5, averagePlaqueCounts([">50", "30"], 53))

    def testgreaterThanNdNumber(self):
        """
        Correct result if one value is nd, one's a > and one's a number.
        """
        self.assertEqual(41.5, averagePlaqueCounts([">50", "30", "nd"], 53))


class TestCalculateFractionInfectivity(TestCase):
    """
    Tests for the calculateFractionInfectivity function.
    """

    def testNd(self):
        """
        If the value is a nd, return nd.
        """
        self.assertEqual("nd", calculateFractionInfectivity("nd", 53))

    def testNeutPercGreaterThan(self):
        """
        If value is >, return correct fraction. This assumes that the >
        means no neutralisation.
        """
        self.assertEqual(
            0, calculateFractionInfectivity(">50", 53, neutralisationPercent=True)
        )

    def testNeutPercNumber(self):
        """
        Calculate correct neutralisation percent if number is given.
        """
        self.assertEqual(
            50, calculateFractionInfectivity("26", 52, neutralisationPercent=True)
        )

    def testFracGreaterThan(self):
        """
        Correct fraction infectivity if > number is given.
        """
        self.assertEqual(
            1, calculateFractionInfectivity(">50", 53, neutralisationPercent=False)
        )

    def testFracNumber(self):
        """
        Correct fraction infectivity with real number.
        """
        self.assertEqual(
            0.5, calculateFractionInfectivity("26", 52, neutralisationPercent=False)
        )

    def testFracNumberOther(self):
        """
        Correct fraction infectivity with real number.
        """
        self.assertEqual(
            0.25, calculateFractionInfectivity("13", 52, neutralisationPercent=False)
        )

    def testFracNumberInt(self):
        """
        Correct fraction infectivity with real number if given as int.
        """
        self.assertEqual(
            0.25, calculateFractionInfectivity(13, 52, neutralisationPercent=False)
        )


class TestConvertRawCountsToNeutcurveDf(TestCase):
    """
    Tests for the convertRawCountsToNeutcurveDf function.
    """

    def testData(self):
        """
        Tests for the convertRawCountsToNeutcurveDf function.
        """
        data = parseTiterExcel(TESTDATA, 10, 3, VIRUSES, addAverage=True)
        d = convertRawCountsToNeutcurveDf(data)

        # Test header
        self.assertEqual(
            ["concentration", "fraction infectivity", "virus", "serum", "replicate"],
            list(d.columns),
        )

        # Test length
        self.assertEqual(4360, d.shape[0])

        # Test specific cases
        self.assertEqual(
            [0.00625, 0.64, "SARS-CoV-2_WT (984)", "1.1", "1"], list(d.loc[0])
        )

        self.assertEqual(
            [0.0015625, 0.0, "SARS-CoV-2_Delta (25853_23)", "4.1", "1"],
            list(d.loc[1036]),
        )

    def testDataCustomTiterSteps(self):
        """
        Test result when custom titersteps are given.
        """
        ts = [
            "1:20",
            "1:40",
            "1:54",
            "1:108",
            "1:216",
            "1:432",
            "1:864",
            "1:1728",
            "1:5120",
        ]

        data = parseTiterExcel(
            TESTDATA, 10, 3, VIRUSES, addAverage=True, serumColumnsNames=ts
        )

        # Assert that column headings are correct
        columns = [
            "sample ID",
            "Infektion Hamster",
            "Antigen",
            "Viruskontrolle",
            "Replicate",
            "1:20",
            "1:40",
            "1:54",
            "1:108",
            "1:216",
            "1:432",
            "1:864",
            "1:1728",
            "1:5120",
        ]
        self.assertEqual(columns, list(data.columns))

        d = convertRawCountsToNeutcurveDf(data, customTiterSteps=ts)

        # Test specific cases
        self.assertEqual(
            [0.009259259259259259, 0.64, "SARS-CoV-2_WT (984)", "1.1", "1"],
            list(d.loc[0]),
        )

        self.assertEqual(
            [0.0023148148148148147, 0.0, "SARS-CoV-2_Delta (25853_23)", "4.1", "1"],
            list(d.loc[1036]),
        )


class TestConvertRawCountsToDiscreteDf(TestCase):
    """
    Tests for the convertRawCountsToDiscreteDf function.
    """

    def testData(self):
        """
        Tests for the convertRawCountsToDiscreteDf function.
        """
        data = parseTiterExcel(TESTDATA, 10, 3, VIRUSES, addAverage=True)
        d = convertRawCountsToDiscreteDf(data)

        # Check correct columns
        self.assertEqual(
            [
                "sample ID",
                "Infektion Hamster",
                "Antigen",
                "1:20",
                "1:40",
                "1:80",
                "1:160",
                "1:320",
                "1:640",
                "1:1280",
                "1:2560",
                "1:5120",
            ],
            list(d.columns),
        )

        # Check that the correct viruses are present in the dataframe
        self.assertEqual(set(VIRUSES), set(data["Antigen"]))

        # Check that the correct sera are present in the dataframe
        self.assertEqual(set(SERA), set(data["sample ID"]))

        # Check some individual cases
        self.assertEqual(
            [
                "1.1",
                "BA.2 Isolate 26729_2",
                "BF.7",
                "nd",
                "nd",
                98.91304347826087,
                88.26086956521739,
                61.95652173913043,
                32.39130434782609,
                0.0,
                0.0,
                "nd",
            ],
            list(d.loc[0]),
        )

        self.assertEqual(
            [
                "1.1",
                "BA.2 Isolate 26729_2",
                "SARS-CoV-2_Delta (25853_23)",
                89.50902224087285,
                68.52706672261854,
                11.875786823331921,
                -11.2043642467478,
                -17.49895090222408,
                -36.382710868652964,
                -36.382710868652964,
                -17.49895090222408,
                -13.302559798573242,
            ],
            list(d.loc[12]),
        )

        self.assertEqual(
            [
                "3.2",
                "Beta (22131)",
                "SARS-CoV-2_Delta (25853_23)",
                "nd",
                "nd",
                "nd",
                97.90180444817457,
                68.52706672261854,
                32.85774234158623,
                11.875786823331921,
                -34.28451531682754,
                -44.7754930759547,
            ],
            list(d.loc[165]),
        )

    def testConvertRawCountsToDiscreteDfCustomTiterSteps(self):
        """
        Result must be correct if custom titer steps are given
        """
        ts = [
            "1:20",
            "1:40",
            "1:54",
            "1:108",
            "1:216",
            "1:432",
            "1:864",
            "1:1728",
            "1:5120",
        ]

        data = parseTiterExcel(
            TESTDATA, 10, 3, VIRUSES, addAverage=True, serumColumnsNames=ts
        )
        d = convertRawCountsToDiscreteDf(data, customTiterSteps=ts)

        # Check some individual cases
        self.assertEqual(
            [
                "1.1",
                "BA.2 Isolate 26729_2",
                "BF.7",
                "nd",
                "nd",
                98.91304347826087,
                88.26086956521739,
                61.95652173913043,
                32.39130434782609,
                0.0,
                0.0,
                "nd",
            ],
            list(d.loc[0]),
        )

        self.assertEqual(
            [
                "1.1",
                "BA.2 Isolate 26729_2",
                "SARS-CoV-2_Delta (25853_23)",
                89.50902224087285,
                68.52706672261854,
                11.875786823331921,
                -11.2043642467478,
                -17.49895090222408,
                -36.382710868652964,
                -36.382710868652964,
                -17.49895090222408,
                -13.302559798573242,
            ],
            list(d.loc[12]),
        )

        self.assertEqual(
            [
                "3.2",
                "Beta (22131)",
                "SARS-CoV-2_Delta (25853_23)",
                "nd",
                "nd",
                "nd",
                97.90180444817457,
                68.52706672261854,
                32.85774234158623,
                11.875786823331921,
                -34.28451531682754,
                -44.7754930759547,
            ],
            list(d.loc[165]),
        )


class TestGetPRNTContinuous(TestCase):
    """
    Tests for the getPRNTContinuous function
    """

    def testFullCurve(self):
        """
        The correct titers must be returned if the titercurve has an
        inflection point.
        """
        d = parseTiterExcel(TESTDATA, 10, 3, VIRUSES, addAverage=True)

        dataSubset = d.loc[
            (d["Antigen"] == "rSARS-CoV-2_B.1+E484K") & (d["sample ID"] == "1.3")
        ]

        dSubset = convertRawCountsToNeutcurveDf(dataSubset)
        self.assertEqual("203.52", getPRNTContinuous(dSubset)[0])
        self.assertEqual("203.21", getPRNTContinuous(dSubset, fixtop=False)[0])
        self.assertEqual("205.18", getPRNTContinuous(dSubset, fixbottom=False)[0])
        self.assertEqual(
            "205.17", getPRNTContinuous(dSubset, fixtop=False, fixbottom=False)[0]
        )

    def testGreaterThan(self):
        """
        The correct titer must be returned if the serum neutralises in all
        dilutions.
        """
        d = parseTiterExcel(TESTDATA, 10, 3, VIRUSES, addAverage=True)

        dataSubset = d.loc[
            (d["Antigen"] == "SARS-CoV-2_Alpha (21528)") & (d["sample ID"] == "2.1")
        ]

        dSubset = convertRawCountsToNeutcurveDf(dataSubset)
        self.assertEqual(">5120", getPRNTContinuous(dSubset)[0])
        self.assertEqual(">5120", getPRNTContinuous(dSubset, fixtop=False)[0])
        self.assertEqual(">5120", getPRNTContinuous(dSubset, fixbottom=False)[0])
        self.assertEqual(
            ">5120", getPRNTContinuous(dSubset, fixtop=False, fixbottom=False)[0]
        )

    def testLessThan(self):
        """
        The correct titer must be returned if the serum never neutralises.
        """
        d = parseTiterExcel(TESTDATA, 10, 3, VIRUSES, addAverage=True)

        dataSubset = d.loc[
            (d["Antigen"] == "SARS-CoV-2_BA.2 (26729_2)") & (d["sample ID"] == "8.2")
        ]

        dSubset = convertRawCountsToNeutcurveDf(dataSubset)
        self.assertEqual("<160", getPRNTContinuous(dSubset)[0])
        self.assertEqual("<160", getPRNTContinuous(dSubset, fixtop=False)[0])
        self.assertEqual("<160", getPRNTContinuous(dSubset, fixbottom=False)[0])
        self.assertEqual(
            "<160", getPRNTContinuous(dSubset, fixtop=False, fixbottom=False)[0]
        )


class TestAveragePlaqueCounts(TestCase):
    """
    Tests for the averagePlaqueCounts function.
    """

    def testAllNds(self):
        """
        Test when all plaque counts are nd.
        """
        counts = ["nd", "nd", "nd"]
        result = averagePlaqueCounts(counts, 50)
        self.assertEqual(result, "nd")

    def testAllGreaterThan(self):
        """
        Test when all plaque counts are >.
        """
        counts = [">50", ">50", ">50"]
        result = averagePlaqueCounts(counts, 50)
        self.assertEqual(result, ">50")

    def testGreaterThanNds(self):
        """
        Test when plaque counts are mixture of nd and >.
        """
        counts = ["nd", "nd", "nd", ">50", ">50", ">50"]
        result = averagePlaqueCounts(counts, 50)
        self.assertEqual(result, ">50")

    def testAveraging(self):
        """
        Test the averaging of other titers.
        """
        counts = ["20", "20", ">50"]
        result = averagePlaqueCounts(counts, 50)
        self.assertEqual(result, 30.0)


class TestCalculateFractionInfectivity(TestCase):
    """
    Tests for the calculateFractionInfectivity function.
    """

    def testNd(self):
        """
        Test when the number of plaques is nd.
        """
        result = calculateFractionInfectivity("nd", 50)
        self.assertEqual(result, "nd")

    def testGreaterThan(self):
        """
        Test when the number of plaques is >.
        """
        result1 = calculateFractionInfectivity(">50", 50)
        self.assertEqual(result1, 1)

        result2 = calculateFractionInfectivity(">50", 50, neutralisationPercent=True)
        self.assertEqual(result2, 0)

    def testNeutralisationPercentTrue(self):
        """
        Test when neutralisationPercent=True.
        """
        result = calculateFractionInfectivity(40, 50, neutralisationPercent=True)
        self.assertEqual(result, 20)

    def testNeutralisationPercentFalse(self):
        """
        Test when neutralisationPercent=Fale.
        """
        result = calculateFractionInfectivity(40, 50)
        self.assertEqual(result, 0.8)
