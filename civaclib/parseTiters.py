import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from warnings import warn
import sys

from math import log2

import neutcurve

from .common import titerSteps, LIMIT


def averagePlaqueCounts(counts, vk):
    """
    Calculate the average number of plaques from a list of plaque counts.
    Adequately deal with 'nd's, '>'thans etc.

    @param counts: a C{list} of plaque counts to average. These are strings.
    @param vk: the viruskontrolle.
    """
    countsSet = set(counts)

    if countsSet == {"nd"}:
        return "nd"
    elif countsSet == {">50"}:
        return ">50"
    elif countsSet == {">50", "nd"}:
        return ">50"
    else:
        toAverage = []
        for c in counts:
            if str(c).startswith(">"):
                # The exact number of plaques couldn't be determined but is
                # higher than the seeding dose.
                toAverage.append(float(vk))
            elif c == "e":
                # 'e' indicates that the number of plaques corresponds to the
                # seeding dose (used when more than two dilutions in a series
                # have the same pattern). Used starting meeting on 2023-11-17.
                toAverage.append(float(vk))
            elif c == "nd":
                # No titration was carried out.
                continue
            else:
                # This is a plaque count.
                toAverage.append(float(c))
        return np.mean(toAverage)


def parseTiterExcel(
    fileName,
    rowTiterStart,
    serumColumnsStart,
    viruses,
    serumColumnsNames=None,
    addAverage=False,
):
    """
    Parse the excel file with titers from Felix and Marie. The file has the
    following columns: 'sample ID', 'Infektion Hamster', a column with
    information about repeats. Then follows a set of columns that's repeated
    for each virus that the sera were titrated against consisting of
    'Viruskontrolle (VK)', and the dilutions. The first rows contain
    information about the Viruskontrolle. The information about the actual
    titers starts on the given by the index in 'rowTiterStart'.

    @param fileName: the C{str} filename of the excel file.
    @param rowTiterStart: the C{int} index of the row on which the titrations
        start.
    @param serumColumnsStart: the C{int} index of the column of the first
        virus.
    @param viruses: a C{list} of virus names that the sera were titrated
        against.
    @param serumColumnsNames: a C{list} of column headers present for all
        viruses.
    @param addAverage: if C{True}: add the average of the replicates per
        serum/ag to the dataframe.

    @return: a C{pandas.DataFrame} with raw plaque counts parsed out.
    """
    serumColumnsNames = serumColumnsNames or titerSteps

    # Read in the filename
    d = pd.read_csv(fileName)

    # Sanity checks
    missingViruses = set(v for v in d.columns if "Unnamed" not in v) - set(viruses)
    if missingViruses:
        warn(
            f"These viruses specified by you are not present in the "
            f"file: {missingViruses}",
            RuntimeWarning,
        )

    # Change the index
    d.columns = d.iloc[0]

    # Get the indices of the titers
    firstTiterFound = 0

    indices = {}

    for i, colname in enumerate(d.columns):
        if colname == "1:20":
            virusName = viruses[firstTiterFound]
            indices[virusName] = [i, i + len(serumColumnsNames) - 1]
            firstTiterFound += 1

    parsedTable = []
    for virusIndex, virus in enumerate(viruses):
        # Get the indices of the first and column that we're considering for
        # this particular virus
        virusColumnInfo = indices[virus]

        # Make a note of the index of the first and last columns with plaque
        # counts
        countStart = virusColumnInfo[0]
        countEnd = virusColumnInfo[1]

        # Make sure the indices are correct
        assert d.columns[countStart] == serumColumnsNames[0]
        assert d.columns[countEnd] == serumColumnsNames[-1]

        # Loop through the rows, and if it's a row related to raw plaque
        # counts, keep them.
        for i, row in d.iterrows():
            if i >= rowTiterStart:
                rowList = list(row)
                if rowList[0] == rowList[0]:
                    # If the first element in the row isn't nan, the first element is
                    # the name of the serum, the second element the name of the virus.
                    serumName = rowList[0]
                    serumVirusName = rowList[1]
                if rowList[2] in ["PFU Ansatz 1", "PFU Ansatz 2"]:
                    # Get the information about the Viruskontrolle
                    averageVirusKontrolle = None

                    # Loop over the possible columns to find the column with the viruskontrolle
                    for offset in (1, 2, 3):
                        if rowList[countStart - offset] == "x":
                            averageVirusKontrolle = float(
                                d.iat[rowTiterStart - 1, countStart - offset]
                            )
                            break
                    # Check if averageVirusKontrolle was found
                    if averageVirusKontrolle is None:
                        raise ValueError("No Viruskontrolle found for virus %s" % virus)

                    # Add the data
                    rowData = [
                        serumName,
                        serumVirusName,
                        virus,
                        averageVirusKontrolle,
                        rowList[2],
                    ] + rowList[countStart : countEnd + 1]
                    parsedTable.append(rowData)

    countData = pd.DataFrame(
        parsedTable,
        columns=[
            "sample ID",
            "Infektion Hamster",
            "Antigen",
            "Viruskontrolle",
            "Replicate",
        ]
        + serumColumnsNames,
    )

    if not addAverage:
        return countData
    else:
        # Add information about the average plaque counts of all replicates
        # per serum/virus pair.
        averages = []
        for group in countData.groupby(by=["sample ID", "Antigen"]):
            groupData = []
            groupData.extend(
                [
                    list(group[1]["sample ID"])[0],
                    list(group[1]["Infektion Hamster"])[0],
                    list(group[1]["Antigen"])[0],
                    list(group[1]["Viruskontrolle"])[0],
                    "Average",
                ]
            )

            for titerStep in serumColumnsNames:
                counts = list(group[1][titerStep])
                avg = averagePlaqueCounts(counts, list(group[1]["Viruskontrolle"])[0])
                groupData.append(avg)

            averages.append(groupData)

        averagesDf = pd.DataFrame(
            averages,
            columns=[
                "sample ID",
                "Infektion Hamster",
                "Antigen",
                "Viruskontrolle",
                "Replicate",
            ]
            + serumColumnsNames,
        )
        return pd.concat([countData, averagesDf], axis=0, ignore_index=True)


def calculateFractionInfectivity(nr, vk, neutralisationPercent=False):
    """
    Given the number of plaques (nr) and the number of plaques in the control
    well, calculate the remaining infectivity after neutralisation.

    @param nr: the number of plaques with sera. These are floats if plaques could be
        counted, or strings starting with '>' if there were too many plaques to be
        counted.
    @param vk: the number of plaques without sera (control).
    @param neutralisationPercent: if C{bool} True, returns the percentage of
        plaques that were neutralised. If C{False} returns the fraction of
        infectivity remaining.
    """
    if nr == "nd":
        return "nd"
    if neutralisationPercent:
        if str(nr).startswith(">"):
            # Assumes that '>' indicates no neutralisation
            return 0
        elif nr == "e":
            # 'e' means that the number of plaques corresponds to the seeding
            # dose and thus no neutralization.
            return 0
        else:
            return 100 - (float(nr) / vk) * 100
    else:
        if str(nr).startswith(">"):
            # Assumes that '>' indicates no neutralisation
            return 1
        elif nr == "e":
            # 'e' means that the number of plaques corresponds to the seeding
            # dose and thus no neutralization.
            return 1
        else:
            return float(nr) / vk


def convertRawCountsToNeutcurveDf(rawDf, dropAverages=True, customTiterSteps=None):
    """
    Convert the dataframe from the raw count data into the format
    required by neutcurve.

    @param rawDf: a C{pandas.DataFrame} from reading the parsed raw count
        data.
    @param dropAverages: whether rows with average plaque counts should be
        dropped.
    @param customTiterSteps: A C{list} of dilutions that were titrated.
    """
    customTiterSteps = customTiterSteps or titerSteps

    concentrations = [1 / int(step[2:]) for step in customTiterSteps]

    colNames = [
        "sample ID",
        "Infektion Hamster",
        "Antigen",
        "Viruskontrolle",
        "Replicate",
    ] + customTiterSteps

    rawDf.columns = colNames

    rawDf = rawDf.copy(deep=True)

    for i, titerStep in enumerate(customTiterSteps):
        rawDf[concentrations[i]] = [
            calculateFractionInfectivity(row[titerStep], row["Viruskontrolle"])
            for i, row in rawDf.iterrows()
        ]

    if dropAverages:
        rawDf = rawDf.loc[rawDf["Replicate"] != "Average"]

    # Reformat the data so it has the right columns for neutcurve.
    concentration = []
    fractionInfectivity = []
    virus = []
    serum = []
    replicate = []

    for i, row in rawDf.iterrows():
        for conc in concentrations:
            if row[conc] != "nd":
                concentration.append(conc)
                fractionInfectivity.append(row[conc])
                virus.append(row["Antigen"])
                serum.append(row["sample ID"])
                replicate.append(row["Replicate"][-1:])

    neutcurveData = pd.DataFrame(
        {
            "concentration": concentration,
            "fraction infectivity": fractionInfectivity,
            "virus": virus,
            "serum": serum,
            "replicate": replicate,
        }
    )

    return neutcurveData


def convertRawCountsToDiscreteDf(rawDf, customTiterSteps=None):
    """
    Convert the dataframe from the raw count data into the format
    required to calculate discrete PRNT50 titers.

    @param rawDf: a C{pandas.DataFrame} from reading the parsed raw count
        data.
    """
    customTiterSteps = customTiterSteps or titerSteps

    # Only look at rows with average plaque counts.
    if rawDf.shape[0] > 1 and "Average" in list(rawDf["Replicate"]):
        rawDf = rawDf.loc[rawDf["Replicate"] == "Average"]
    elif rawDf.shape[0] > 1 and "Average" not in list(rawDf["Replicate"]):
        raise ('Dataframe must contain "Average" replicate or only one replicate.')

    reformatted = []
    for i, row in rawDf.iterrows():
        rowInfo = []
        rowInfo.extend([row["sample ID"], row["Infektion Hamster"], row["Antigen"]])

        for titerStep in customTiterSteps:
            rowInfo.append(
                calculateFractionInfectivity(
                    row[titerStep], row["Viruskontrolle"], neutralisationPercent=True
                )
            )
        reformatted.append(rowInfo)

    columnNames = ["sample ID", "Infektion Hamster", "Antigen"] + customTiterSteps

    return pd.DataFrame(reformatted, columns=columnNames)


def getPRNTDiscrete(rawCounts, limit=None, steps=None, nd="<20"):
    """
    From a dictionary of titer steps and percent reduction of number of
    plaques, get the actual titer.

    @param rawCounts: a C{dict} where the keys are dilution steps ('1:20')
        and the values are the percent reduction in number of plaques.
    @param limit: the C{int} level of sensitivity. Usually 50 or 90.
    @param steps: a C{list} of dilution levels that were tested. This should
        be ordered from lowest to highest dilution.
    @param nd: the lowest titer level. Usually '<20'.
    """
    limit = limit or LIMIT
    steps = steps or titerSteps

    # If no titrations have been done, return *
    if set(rawCounts.values()) == {"nd"}:
        return "*"

    assert list(rawCounts) == steps, "Data dilutions do not match given dilutions."

    prnt = nd
    titerFound = False
    first = True

    for i, step in enumerate(steps[::-1]):
        titerAtDilution = rawCounts[step]
        if titerAtDilution != "nd":
            if float(titerAtDilution) > limit:
                if first:
                    titerFound = True
                    # The first titer (highest dilution) in the dilution
                    # series is already below the limit. Set to > than.
                    if float(titerAtDilution - 5) > limit:
                        prnt = ">" + steps[::-1][i][2:]
                    else:
                        prnt = steps[::-1][i][2:]
                else:
                    titerAtPreviousDilution = rawCounts[steps[::-1][i - 1]]
                    if abs(titerAtPreviousDilution - limit) < abs(
                        titerAtDilution - limit
                    ):
                        prnt = steps[::-1][i - 1][2:]
                    else:
                        prnt = steps[::-1][i][2:]
                    titerFound = True
                    first = False
                break
            else:
                first = False

    if not titerFound:
        # This means that the number of plaques never dropped below the level
        # of sensitivity. The titer should be expressed as greater than the
        # highest dilution tested.
        # Therefore, step through the dilution steps from lowest to highest
        # looking for the first time a non-nd titer is encountered.
        for step in steps:
            titerAtDilution = rawCounts[step]
            if titerAtDilution != "nd":
                if float(titerAtDilution + 5) < limit:
                    # Use a cut-off to prevent titers that are close to the
                    # cut-off to be nd
                    prnt = "<" + step[2:]
                    break
                elif float(titerAtDilution) < limit:
                    prnt = step[2:]
                    break
                elif float(titerAtDilution) > limit:
                    print(f"Something went wrong for step {step}.", file=sys.stderr)
                    sys.exit(1)

    return prnt


def getPRNTContinuous(data, limit=50, fixtop=True, fixbottom=True, interpolate=False):
    """
    Get continuous PRNT titers. This uses the neutcurve package by the Bloom
    lab.

    @param limit: the C{int} level of sensitivity. Usually 0.5 or 0.9.
    @param fixtop: Fix the top of the neutralisation curve at 1.
    @param fixbottom: Fix the bottom of the neutralisation curve at 0.
    @param interpolate: If C{True} interpolate titers that are out of bounds of
        the dilutions tested. IF USING THIS, MAKE SURE TO CHECK THE
        INTERPOLATED TITER AGAINST THE NEUTRALISATION CURVE TO MAKE SURE THE
        TITER STILL MAKES SENSE!!
    """
    limit = limit / 100

    # If no titrations have ben done, return *
    if data.shape[0] == 0:
        return "*", "*"

    if fixtop and fixbottom:
        curve = neutcurve.HillCurve(
            data["concentration"], data["fraction infectivity"], fitlogc=False
        )

    elif fixtop and fixbottom is False:
        curve = neutcurve.HillCurve(
            data["concentration"],
            data["fraction infectivity"],
            fitlogc=False,
            fixbottom=False,
        )

    elif fixtop is False and fixbottom:
        curve = neutcurve.HillCurve(
            data["concentration"],
            data["fraction infectivity"],
            fitlogc=False,
            fixtop=False,
        )

    elif fixtop is False and fixbottom is False:
        curve = neutcurve.HillCurve(
            data["concentration"],
            data["fraction infectivity"],
            fitlogc=False,
            fixtop=False,
            fixbottom=False,
        )

    if curve.icXX(limit):
        # The titer is within one of the performed dilutions of the assay.
        prnt = f"{1/curve.icXX(limit):.2f}"
    else:
        # The titer is non-detectable, either as > or <.
        if curve.icXX_str(limit).startswith("<"):
            cutoffs = {
                "<0.000195": {
                    "limit": 10240,
                    "newTiter": ">5120",
                },
                "<0.000391": {
                    "limit": 5120,
                    "newTiter": ">2560",
                },
                "<0.000781": {
                    "limit": 2560,
                    "newTiter": ">1280",
                },
                "<0.00156": {
                    "limit": 1280,
                    "newTiter": ">640",
                },
                "<0.00313": {
                    "limit": 640,
                    "newTiter": ">320",
                },
                "<0.00625": {
                    "limit": 320,
                    "newTiter": ">160",
                },
            }
            # This is a > non-detectable titer.
            # The neutcurve program doesn't give enough decimal points to get to
            # the correct discrete titer value. Therefore do it by hand.
            neutcurveLimit = curve.icXX_str(limit)
            if neutcurveLimit in cutoffs:
                if interpolate:
                    prnt = interpolateContinuous(curve, limit)
                    if prnt > cutoffs[neutcurveLimit]["limit"]:
                        # I don't want to interpolate too far out of range.
                        prnt = cutoffs[neutcurveLimit]["newTiter"]
                    else:
                        prnt = f"{1/prnt:.2f}"
                else:
                    prnt = cutoffs[neutcurveLimit]["newTiter"]
            else:
                # This is a <?? non-detectable titer, but the upper bound isn't
                # in the list above.
                print("unknown upper bound", curve.icXX_str(limit))
                prnt = f">{int(1/float(curve.icXX_str(limit)[1:]))}"
        elif curve.icXX_str(limit).startswith(">"):
            # This is a < non-detectable titer.
            prnt = f"<{int(1/float(curve.icXX_str(limit)[1:]))}"

    return prnt, curve


def interpolateContinuous(curve, limit):
    """
    Interpolate the curve returned by `neutcurve.HillCurve` outside the
    tested dilutions. ONLY USE THIS WITH CAUTION AND CHECK RESULTS AGAINST
    THE TITER CURVE TO MAKE SURE THE TITER MAKES SENSE!!!
    This code is slightly adapted from the neutcurve.HillCurve.icXX
    function.

    @param curve: the C{neutcurve.HillCurve} curve fit.
    @param limit: the C{int} level of sensitivity. Usually 0.5 or 0.9.
    """
    fracinf = 1 - limit

    icXX = curve.midpoint * ((curve.top - fracinf) / (fracinf - curve.bottom)) ** (
        1.0 / curve.slope
    )

    if not (curve.cs[0] <= icXX <= curve.cs[-1]):
        # Return the interpolated titer if the titer is out of bounds.
        return icXX
    else:
        # I could just return the interpolated titer here anyway. But I want
        # to make sure this function is only used when absolutely necessary.
        raise ValueError(
            "You shouldn't be calling 'interpolateContinuous' "
            "with a titer that falls into the tested dilutions."
        )


def getAllTiters(
    data, plot=False, ax=False, interpolate=False, limit=50, customTiterSteps=None
):
    """
    Return discrete and continuous titers. If requested, plot the titer curve.

    @param data: the C{pandas.dataframe} returned by {parseTiterExcel} and
        subsetted to the viruses and sera to get the titers for/plot.
    @param plot: if C{True} generate a plot of the titer curves.
    @param ax: if C{True} an axes object to use for plotting.
    @param interpolate: If C{True} interpolate titers that are out of bounds of
        the dilutions tested. IF USING THIS, MAKE SURE TO CHECK THE
        INTERPOLATED TITER AGAINST THE NEUTRALISATION CURVE TO MAKE SURE THE
        TITER STILL MAKES SENSE!!
    """
    customTiterSteps = customTiterSteps or titerSteps
    titerStepsNumbers = [int(step[2:]) for step in customTiterSteps]
    finalDilution = titerStepsNumbers[-1] * 2

    if ax and not plot:
        raise "Plotting not specified but axes object is given."

    continuousData = convertRawCountsToNeutcurveDf(
        data, customTiterSteps=customTiterSteps
    )
    discreteData = convertRawCountsToDiscreteDf(data, customTiterSteps=customTiterSteps)
    if continuousData.shape[0] == 0:
        return ("*", "*", "*", "*", "*")
    assert discreteData.shape[0] == 1, "Subset the dataframe given"

    prnt50ContFixtopFixbottom, curveFixtopFixbottom = getPRNTContinuous(
        continuousData,
        fixtop=True,
        fixbottom=True,
        interpolate=interpolate,
        limit=limit,
    )

    prnt50ContFixtop, curveFixTop = getPRNTContinuous(
        continuousData,
        fixtop=True,
        fixbottom=False,
        interpolate=interpolate,
        limit=limit,
    )

    prnt50ContFixbottom, curveFixbottom = getPRNTContinuous(
        continuousData,
        fixtop=False,
        fixbottom=True,
        interpolate=interpolate,
        limit=limit,
    )

    prnt50Cont, curve = getPRNTContinuous(
        continuousData,
        fixtop=False,
        fixbottom=False,
        interpolate=interpolate,
        limit=limit,
    )

    for i, row in discreteData.iterrows():
        d = {titerStep: row[titerStep] for titerStep in customTiterSteps}

        prnt50discrete = getPRNTDiscrete(
            d, limit=limit, steps=customTiterSteps, nd="<20"
        )

    if not plot:
        return (
            prnt50discrete,
            prnt50ContFixtopFixbottom,
            prnt50ContFixtop,
            prnt50ContFixbottom,
            prnt50Cont,
        )
    else:
        if not ax:
            fig, ax = plt.subplots()

        continuousDataRep1 = continuousData.loc[(continuousData["replicate"] == "1")]
        continuousDataRep2 = continuousData.loc[(continuousData["replicate"] == "2")]

        ax.plot(
            continuousDataRep1["concentration"],
            continuousDataRep1["fraction infectivity"],
            "o",
            alpha=0.5,
        )
        ax.plot(
            continuousDataRep2["concentration"],
            continuousDataRep2["fraction infectivity"],
            "s",
            alpha=0.5,
        )

        fittedFixtopFixbottom = curveFixtopFixbottom.dataframe()
        ax.plot(
            fittedFixtopFixbottom["concentration"],
            fittedFixtopFixbottom["fit"],
            "-",
            color="black",
        )

        fittedFixTop = curveFixTop.dataframe()
        ax.plot(
            fittedFixTop["concentration"],
            fittedFixTop["fit"],
            "-",
            color="red",
            linewidth=0.5,
        )

        fittedFixbottom = curveFixbottom.dataframe()
        ax.plot(
            fittedFixbottom["concentration"],
            fittedFixbottom["fit"],
            "-",
            color="blue",
            linewidth=0.5,
        )

        fitted = curve.dataframe()
        ax.plot(
            fitted["concentration"], fitted["fit"], "-", color="grey", linewidth=0.5
        )

        for y in (0.5, 0.1, 0.9):
            ax.hlines(
                y=y,
                xmin=1 / 10,
                xmax=1 / 81920,
                color="grey",
                linestyle="--",
                linewidth=0.5,
            )
        for y in (0, 1):
            ax.hlines(
                y=y,
                xmin=1 / 10,
                xmax=1 / 81920,
                color="grey",
                linestyle="-",
                linewidth=0.1,
            )

        ax.set_xscale("log")

        srNames = {
            "1.1": "1.1 (BA.2-2)",
            "1.2": "1.2 (BA.2-2)",
            "1.3": "1.3 (BA.2-2)",
            "2.1": "2.1 (Alpha)",
            "2.2": "2.1 (Alpha)",
            "3.1": "3.1 (Beta)",
            "3.2": "3.2 (Beta)",
            "3.3": "3.3 (Beta)",
            "4.1": "4.1 (Delta)",
            "4.2": "4.2 (Delta)",
            "4.3": "4.3 (Delta)",
            "5.1": "5.1 (BA.1)",
            "5.2": "5.2 (BA.1)",
            "5.3": "5.3 (BA.1)",
            "6.1": "6.1 (BA.2-12)",
            "6.2": "6.2 (BA.2-12)",
            "6.3": "6.3 (BA.2-12)",
            "7.1": "7.1 (B.1 + E484K)",
            "7.2": "7.2 (B.1 + E484K)",
            "7.3": "7.3 (B.1 + E484K)",
            "8.1": "8.1 (WT)",
            "8.2": "8.2 (WT)",
            "8.3": "8.3 (WT)",
            "9.1": "9.1 (BA.5)",
            "9.2": "9.2 (BA.5)",
            "9.3": "9.3 (BA.5)",
            "12SE0030": "12SE0030 (XBB.2)",
            "12SE0031": "12SE0031 (XBB.2)",
            "12SE0032": "12SE0032 (XBB.2)",
        }

        ax.set_title(
            f'{srNames[list(data["sample ID"])[0]]} vs '
            f'{list(data["Antigen"])[0]}\nPRNT{limit}: disc: '
            f"{prnt50discrete}, contFixBoth: "
            f"{prnt50ContFixtopFixbottom},\n"
            f"contFixTop: {prnt50ContFixtop}, contFixBottom: "
            f"{prnt50ContFixbottom}, cont: {prnt50Cont}"
        )

        ax.set_xticks([])
        ax.set_xticks([1 / tsn for tsn in titerStepsNumbers] + [1 / finalDilution])

        ax.set_xticklabels(customTiterSteps + [f"1:{finalDilution}"], rotation=90)
        ax.set_ylim(-0.5, 1.5)
        ax.set_xlim(1 / 5120, 1 / 20)

        return (
            prnt50discrete,
            prnt50ContFixtopFixbottom,
            prnt50ContFixtop,
            prnt50ContFixbottom,
            prnt50Cont,
        )


def plotTiterCurves(rawCounts, steps=None, ax=None, limit=50, title="", color="black"):
    steps = steps or titerSteps

    assert list(rawCounts) == steps, "Data dilutions do not match given dilutions."

    prnt50 = getPRNTDiscrete(rawCounts, limit=limit, steps=steps, nd="<160")

    if not ax:
        fig, ax = plt.subplots()

    titers = []
    levelIndex = []

    for i, step in enumerate(steps):
        titerAtStep = rawCounts[step]
        titers.append(float(titerAtStep) if titerAtStep != "nd" else np.nan)
        levelIndex.append(i)

    ax.plot(levelIndex, titers, "-o", color=color)
    ax.set_xticks(levelIndex)
    ax.set_xticklabels(steps, rotation=90)
    for y in (0.5, 0.1, 0.9):
        ax.hlines(y=y, xmin=-0.5, xmax=len(levelIndex), color="grey", linestyle="--")
    ax.hlines(
        y=0,
        xmin=-0.5,
        xmax=len(levelIndex),
        color="black",
        linestyle="-",
        linewidth=0.2,
    )
    ax.set_xlim(-0.5, len(levelIndex) - 0.5)
    ax.set_title(f"{title} titer: {prnt50}")
    ax.set_ylim(-20, 100)
    ax.set_ylabel("% reduction in plaques")
    plt.gca().invert_xaxis()
