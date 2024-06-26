# Antigenic cartography using variant-specific hamster sera reveals substantial antigenic variation among Omicron subvariants

This repo contains data and code for the analyses presented in `Antigenic cartography using variant-specific hamster sera reveals substantial antigenic variation among Omicron subvariants`, PNAS, 2024.

## Directory structure

`Rcivaclib` and `civaclib`: R and Python code, respectively.

`bin`: Bin scripts

`data`: Raw data and code to generate antigenic maps.

`figures`: Code to generate the figures shown in the paper, as well as the resulting figures.

`test`: Tests.

## Dependencies

All dependencies for R are listed in the `DEPENDENCIES` file. Apart from the dependencies available through CRAN, the following packages specific to antigenic cartography and titer analysis are required:

##### Racmacs
Available from [github.com/shwilks/Racmacs](github.com/shwilks/Racmacs), or through CRAN.

##### titertools
Available from [github.com/shwilks/titertools](github.com/shwilks/titertools).

Furthermore, for Python, the following libraries need to be installed: `neutcurve`, `pandas`, `numpy`.
