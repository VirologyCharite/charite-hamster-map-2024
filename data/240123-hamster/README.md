# This directory contains the following files

## Raw data
`Hamster_Wdh_PRNT.xlsx` --> Raw data with the titrations of the EG.5.1, BA.2.86, and JN.1 variants, and the XBB.2 sera added.

`PRNT_Hamster_detailliert.csv` --> The raw data used for the analyses. This is the data as of 230420 with the additional titrations of the EG.5.1 and JN.1 variants, and the XBB.2 sera added. For JN.1, the re-titrations from 240123 were used with the exception of serum 1.3 for which the re-titration failed, and the titration from 240105 was used.

`PRNT_Hamster_detailliert-retitration.csv` --> same as `PRNT_Hamster_detailliert.csv`, but with the addition of the original JN.1 titrations from 240105.



## Titertables
`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method discrete --limit 50 --fixedTiterFile data/240123-hamster/adaptations-discrete-50.csv --adaptTiterSteps 230219-xbb2-bn131 > data/240123-hamster/titers-discrete-50.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method discrete --limit 75 --fixedTiterFile data/240123-hamster/adaptations-discrete-75.csv --adaptTiterSteps 230219-xbb2-bn131 > data/240123-hamster/titers-discrete-75.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method discrete --limit 90 --fixedTiterFile data/240123-hamster/adaptations-discrete-90.csv --adaptTiterSteps 230219-xbb2-bn131 > data/240123-hamster/titers-discrete-90.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method discrete --limit 99 --fixedTiterFile data/240123-hamster/adaptations-discrete-99.csv --adaptTiterSteps 230219-xbb2-bn131 > data/240123-hamster/titers-discrete-99.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixtop-fixbottom --limit 50 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixtop-fixbottom-50.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixtop-fixbottom-50.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixtop-fixbottom --limit 75 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixtop-fixbottom-75.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixtop-fixbottom-75.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixtop-fixbottom --limit 90 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixtop-fixbottom-90.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixtop-fixbottom-90.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixtop-fixbottom --limit 99 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixtop-fixbottom-99.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixtop-fixbottom-99.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous --limit 90 --fixedTiterFile data/240123-hamster/adaptations-continuous-90.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-90.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixtop --limit 90 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixtop-90.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixtop-90.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixbottom --limit 90 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixbottom-90.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixbottom-90.csv`

`$ python bin/make-titer-table.py --rawTiters data/240123-hamster/PRNT_Hamster_detailliert.csv --method continuous-fixbottom --limit 90 --fixedTiterFile data/240123-hamster/adaptations-continuous-fixbottom-90-corrected.csv --adaptTiterSteps 230219-xbb2-bn131 --interpolate > data/240123-hamster/titers-continuous-fixbottom-90-corrected.csv`


## Adapted titer files

These files contain adapted titers for titrations where the titer cannot get inferred correctly.

```
adaptations-continuous-90.csv
adaptations-continuous-fixbottom-90-corrected.csv
adaptations-continuous-fixbottom-90.csv
adaptations-continuous-fixtop-90.csv
adaptations-continuous-fixtop-fixbottom-50.csv
adaptations-continuous-fixtop-fixbottom-75.csv
adaptations-continuous-fixtop-fixbottom-90.csv
adaptations-continuous-fixtop-fixbottom-99.csv
adaptations-discrete-50.csv
adaptations-discrete-75.csv
adaptations-discrete-90.csv
adaptations-discrete-99.csv
```

## Maps

Each of the titer csv files have a corresponding map file ending in `*.ace`.
