from os.path import dirname
from pathlib import Path

import civaclib


TOPDIR = dirname(dirname(civaclib.__file__))
TESTDATA = Path(TOPDIR) / "data" / "240123-hamster" / "PRNT_Hamster_detailliert.csv"

titerSteps = [
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

LIMIT = 50


VIRUSES = (
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
)


SERA = [
    "1.1",
    "1.2",
    "1.3",
    "2.1",
    "2.2",
    "3.1",
    "3.2",
    "3.3",
    "4.1",
    "4.2",
    "4.3",
    "5.1",
    "5.2",
    "5.3",
    "6.1",
    "6.2",
    "6.3",
    "7.1",
    "7.2",
    "7.3",
    "8.1",
    "8.2",
    "8.3",
    "9.1",
    "9.2",
    "9.3",
    "12SE0030",
    "12SE0032",
    "12SE0031",
]
