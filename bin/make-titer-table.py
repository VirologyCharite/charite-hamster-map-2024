#!/usr/bin/env python

import sys
import pandas as pd
from warnings import warn

from civaclib.parseTiters import (
    getPRNTDiscrete,
    parseTiterExcel,
    convertRawCountsToDiscreteDf,
    convertRawCountsToNeutcurveDf,
    getPRNTContinuous,
)
from civaclib.common import titerSteps


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Pull out PRNT titers and make a titer table."
    )

    parser.add_argument(
        "--rawTiters",
        help=(
            "The name of the file containing the raw data. See "
            "`data/220706-hamster/PRNT_Hamster_compressed.csv` for a "
            "template of the file format."
        ),
    )

    parser.add_argument(
        "--rowTiterStart",
        default=10,
        type=int,
        help="The index of the row on which the titrations start.",
    )

    parser.add_argument(
        "--serumColumnsStart",
        default=3,
        type=int,
        help="The index of the column of the first virus.",
    )

    parser.add_argument(
        "--limit",
        default=50,
        type=int,
        help="The limit for which titers should be read out.",
    )

    parser.add_argument(
        "--nd",
        default="<20",
        help="The lower titer at below which titers are non-detectable.",
    )

    parser.add_argument(
        "--method",
        default="continuous-fixtop-fixbottom",
        choices=(
            "discrete",
            "continuous-fixtop-fixbottom",
            "continuous-fixtop",
            "continuous-fixbottom",
            "continuous",
        ),
        help="The method used to calculate the titers.",
    )

    parser.add_argument(
        "--interpolate",
        action="store_true",
        help=(
            "If given, interpolate titers that are out of bounds of "
            "the dilutions tested. IF USING THIS, MAKE SURE TO CHECK THE "
            "INTERPOLATED TITERS AGAINST THE NEUTRALISATION CURVE TO MAKE "
            "SURE THE TITERS STILL MAKE SENSE!! You probably also need to "
            "specify a --fixedTiterFile to adjust those titers that are not "
            "interpolated correctly."
        ),
    )

    parser.add_argument(
        "--fixedTiterFile", default=False, help="A csv file of fixed titers."
    )

    parser.add_argument(
        "--adaptTiterSteps",
        default=False,
        help=(
            "Specify if titersteps other than the default titersteps "
            "specified in common.titerSteps should be used. This can be "
            "done by giving keywords corresponding for different "
            "alternative titersteps for particular variant and sera "
            "combinations. For the adaptations to the BN.1.3.1 and XBB.2 "
            'use "230219-xbb2-bn131".'
        ),
    )

    args = parser.parse_args()

    assert args.adaptTiterSteps in {
        "230219-xbb2-bn131",
        False,
    }, "Specify a valid adaptTiterSteps argument."

    if args.fixedTiterFile:
        fixedTiters = pd.read_csv(args.fixedTiterFile)
    else:
        fixedTiters = False

    viruses = [
        "SARS-CoV-2_WT (984)",
        "SARS-CoV-2_Alpha (21528)",
        "SARS-CoV-2_Beta (22131)",
        "SARS-CoV-2_Delta (25853_23)",
        "SARS-CoV-2_BA.1 (26335)",
        "SARS-CoV-2_BA.2 (26729_12)",
        "SARS-CoV-2_BA.2 (26729_2)",
        "rSARS-CoV-2_B.1+E484K",
        "SARS-CoV-2_BA.4",
        "SARS-CoV-2_BA.5",
        "SARS-CoV-2_Mu",
        "BQ.1.18",
        "BF.7",
        "BN.1.3.1",
        "XBB.2",
        "EG.5.1 (V140)",
        "JN.1 (V148)",
    ]
    raw = parseTiterExcel(
        args.rawTiters,
        args.rowTiterStart,
        args.serumColumnsStart,
        viruses,
        addAverage=True,
    )

    titers = []

    if args.method == "discrete":
        if args.interpolate:
            raise ValueError("method=discrete cannot interpolate titers.")
    else:
        if args.interpolate and not args.fixedTiterFile:
            warn(
                "You should probably specify a fixed titer file when "
                "interpolating titers."
            )

    for virus in viruses:
        for serum in set(raw["sample ID"]):
            # Figure out which titersteps to use
            if args.adaptTiterSteps == "230219-xbb2-bn131":
                if serum in {
                    "1.1",
                    "1.2",
                    "1.3",
                    "5.1",
                    "5.2",
                    "5.3",
                    "9.1",
                    "9.2",
                    "9.3",
                } and virus in {"BN.1.3.1", "XBB.2"}:
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
                elif serum in {"8.1", "8.2", "8.3"} and virus == "XBB.2":
                    ts = [
                        "1:20",
                        "1:32",
                        "1:65",
                        "1:130",
                        "1:259",
                        "1:518",
                        "1:1037",
                        "1:2560",
                        "1:5120",
                    ]
                elif serum == "8.3" and virus == "BN.1.3.1":
                    ts = [
                        "1:20",
                        "1:32",
                        "1:65",
                        "1:130",
                        "1:259",
                        "1:518",
                        "1:1037",
                        "1:2560",
                        "1:5120",
                    ]
                else:
                    ts = titerSteps
            else:
                ts = titerSteps

            rawSubset = raw.loc[(raw["Antigen"] == virus) & (raw["sample ID"] == serum)]

            rawSubset.columns = [
                "sample ID",
                "Infektion Hamster",
                "Antigen",
                "Viruskontrolle",
                "Replicate",
            ] + ts

            if args.method == "discrete":
                reformattedDf = convertRawCountsToDiscreteDf(
                    rawSubset, customTiterSteps=ts
                )

                d = {titerStep: list(reformattedDf[titerStep])[0] for titerStep in ts}

                prnt = getPRNTDiscrete(d, limit=args.limit, steps=ts, nd=args.nd)

                titers.append([serum, virus, prnt])

            else:
                # Continuous titers

                reformattedDf = convertRawCountsToNeutcurveDf(
                    rawSubset, customTiterSteps=ts
                )

                if args.method == "continuous-fixtop-fixbottom":
                    fixtop = True
                    fixbottom = True
                elif args.method == "continuous-fixtop":
                    fixtop = True
                    fixbottom = False
                elif args.method == "continuous-fixbottom":
                    fixtop = False
                    fixbottom = True
                else:
                    fixtop = False
                    fixbottom = False

                prnt, curve = getPRNTContinuous(
                    reformattedDf,
                    fixtop=fixtop,
                    fixbottom=fixbottom,
                    interpolate=args.interpolate,
                    limit=args.limit,
                )

                titers.append((serum, virus, prnt))

    titersLong = pd.DataFrame(titers, columns=["serum", "virus", "prnt"])

    # Replace the previously specified fixed titers
    if fixedTiters is not False:
        for i, row in fixedTiters.iterrows():
            origTiter = str(
                titersLong.loc[
                    (titersLong.serum == str(row["serum"]))
                    & (titersLong.virus == row["virus"]),
                    "prnt",
                ].values[0]
            )
            try:
                expectedOrig = f'{row["prntOrig"]:.2f}'
            except ValueError:
                # Deal with '<' and '>' expected original titers
                expectedOrig = row["prntOrig"]

            if origTiter == expectedOrig:
                titersLong.loc[
                    (titersLong.serum == str(row["serum"]))
                    & (titersLong.virus == row["virus"]),
                    "prnt",
                ] = row["prnt"]
            else:
                raise ValueError(
                    f'For {row["serum"]} vs {row["virus"]}, the original '
                    f'titers differ ({origTiter} vs {row["prntOrig"]}).'
                )

    # Convert from long to wide format
    titersWide = pd.pivot(titersLong, index="virus", columns="serum", values="prnt")

    titersWide.to_csv(sys.stdout)
