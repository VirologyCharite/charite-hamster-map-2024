import sys
from os.path import dirname, join

import numpy as np
import pandas as pd
import math

import matplotlib.pyplot as plt

import civaclib

basePath = dirname(dirname(civaclib.__file__))

from civaclib.parseTiters import parseTiterExcel, getAllTiters

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
    "BA.2.86 (V139)",
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

raw = parseTiterExcel(
    join(basePath, "data/240123-hamster/PRNT_Hamster_detailliert.csv"),
    10,
    3,
    viruses,
    addAverage=True,
)

fig, ax = plt.subplots(nrows=29, ncols=17, figsize=(80, 120))

for rowIndex, serum in enumerate(sera):
    for colIndex, virus in enumerate(
        [virus for virus in viruses if virus != "BA.2.86 (V139)"]
    ):

        rawSubset = raw.loc[
            (raw["Antigen"] == virus) & (raw["sample ID"] == str(serum))
        ]

        (
            prnt50discrete,
            prnt50ContFixtopFixbottom,
            prnt50ContFixtop,
            prnt50ContFixbottom,
            prnt50Cont,
        ) = getAllTiters(
            rawSubset, plot=True, ax=ax[rowIndex][colIndex], interpolate=True, limit=90
        )

plt.tight_layout()
plt.savefig(join(basePath, "figures/titercurves/titercurves_90_240123.png"))
