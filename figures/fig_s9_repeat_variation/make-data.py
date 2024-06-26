import numpy as np
import pandas as pd

import math

import sys
from os.path import dirname, join

import civaclib
from civaclib.parseTiters import parseTiterExcel, getAllTiters

# Look at repeat variation between runs

basePath = dirname(dirname(civaclib.__file__))


def convertTiter(titer):
    titer = str(titer)
    if titer == "*":
        titertype = 0
        return np.nan, titertype
    if titer.startswith("<"):
        titertype = 2
        logtiter = np.log2(int(titer[1:]) / 10) - 1
    elif titer.startswith(">"):
        titertype = 2
        logtiter = np.log2(int(titer[1:]) / 10) + 1
    else:
        titertype = 1
        logtiter = np.log2(float(titer) / 10)
    return logtiter, titertype


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

sera = [
    1.1,
    1.2,
    1.3,
    2.1,
    2.2,
    3.1,
    3.2,
    3.3,
    4.1,
    4.2,
    4.3,
    5.1,
    5.2,
    5.3,
    6.1,
    6.2,
    6.3,
    7.1,
    7.2,
    7.3,
    8.1,
    8.2,
    8.3,
    9.1,
    9.2,
    9.3,
    "12SE0030",
    "12SE0031",
    "12SE0032",
]

fixedTiters = pd.read_csv(
    join(basePath, "data/240123-hamster/adaptations-discrete-90.csv")
)
fixedTitersCont = pd.read_csv(
    join(
        basePath,
        "data/240123-hamster/adaptations-continuous-fixbottom-90-corrected.csv",
    )
)

raw = parseTiterExcel(
    join(basePath, "data/240123-hamster/PRNT_Hamster_detailliert.csv"),
    10,
    3,
    viruses,
    addAverage=True,
)

differences = []
for serum in sera:
    for virus in viruses:

        # exclude serum and virus pairs where there are instances where only
        # one repeat was done / countable.
        if f"{virus} {serum}" != "SARS-CoV-2_BA.2 (26729_2) 7.1":

            virusDiffs = [virus, serum]

            rawSubset1 = raw.loc[
                (raw["Antigen"] == virus)
                & (raw["sample ID"] == str(serum))
                & (raw["Replicate"] == "PFU Ansatz 1")
            ]

            (
                prnt50discrete1,
                prnt50ContFixtopFixbottom1,
                prnt50ContFixtop1,
                prnt50ContFixbottom1,
                prnt50Cont1,
            ) = getAllTiters(rawSubset1, plot=False, interpolate=True, limit=90)

            rawSubset2 = raw.loc[
                (raw["Antigen"] == virus)
                & (raw["sample ID"] == str(serum))
                & (raw["Replicate"] == "PFU Ansatz 2")
            ]

            (
                prnt50discrete2,
                prnt50ContFixtopFixbottom2,
                prnt50ContFixtop2,
                prnt50ContFixbottom2,
                prnt50Cont2,
            ) = getAllTiters(rawSubset2, plot=False, interpolate=True, limit=90)

            # Adapt titers
            if len(
                fixedTiters.loc[
                    (fixedTiters["serum"] == serum) & (fixedTiters["virus"] == virus)
                ]["prnt"].values
            ):
                prnt50discrete1 = fixedTiters.loc[
                    (fixedTiters["serum"] == serum) & (fixedTiters["virus"] == virus)
                ]["prnt"].values[0]
                prnt50discrete2 = fixedTiters.loc[
                    (fixedTiters["serum"] == serum) & (fixedTiters["virus"] == virus)
                ]["prnt"].values[0]

            if len(
                fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values
            ):
                prnt50ContFixtopFixbottom1 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]
                prnt50ContFixtop1 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]
                prnt50ContFixbottom1 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]
                prnt50Cont1 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]

                prnt50ContFixtopFixbottom2 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]
                prnt50ContFixtop2 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]
                prnt50ContFixbottom2 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]
                prnt50Cont2 = fixedTitersCont.loc[
                    (fixedTitersCont["serum"] == serum)
                    & (fixedTitersCont["virus"] == virus)
                ]["prnt"].values[0]

            virusDiffs.append(
                convertTiter(prnt50discrete1)[0] - convertTiter(prnt50discrete2)[0]
            )
            virusDiffs.append(
                convertTiter(prnt50ContFixtopFixbottom1)[0]
                - convertTiter(prnt50ContFixtopFixbottom2)[0]
            )
            virusDiffs.append(
                convertTiter(prnt50ContFixtop1)[0] - convertTiter(prnt50ContFixtop2)[0]
            )
            virusDiffs.append(
                convertTiter(prnt50ContFixbottom1)[0]
                - convertTiter(prnt50ContFixbottom2)[0]
            )
            virusDiffs.append(
                convertTiter(prnt50Cont1)[0] - convertTiter(prnt50Cont2)[0]
            )

            differences.append(virusDiffs)

repeatVar = pd.DataFrame(
    differences,
    columns=[
        "Antigen",
        "sample ID",
        "prnt50discrete",
        "prnt50ContFixtopFixbottom",
        "prnt50ContFixtop",
        "prnt50ContFixbottom",
        "prnt50Cont",
    ],
)

repeatVar.to_csv(join(basePath, "figures/fig_s9_repeat_variation/data_240123.csv"))
