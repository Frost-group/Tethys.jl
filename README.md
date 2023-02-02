# Tethys <img src="img/Tethys.svg" width=300 align=right>

[![Build Status](https://travis-ci.org/jarvist/Tethys.jl.svg?branch=master)](https://travis-ci.org/jarvist/Tethys.jl)
[![Coverage Status](https://coveralls.io/repos/jarvist/Tethys.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jarvist/Tethys.jl?branch=master)
[![codecov.io](http://codecov.io/github/jarvist/Tethys.jl/coverage.svg?branch=master)](http://codecov.io/github/jarvist/Tethys.jl?branch=master)
[![docs-latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://jarvist.github.io/Tethys.jl/)

Diagrammatic Monte-Carlo (DiagMC) calculation of Fröhlich polarons. A WIP.

Starting with §3 in [1] - a minimalist bare expansion for the Green function of the Frohlich polaron.

## Inputs
- Electron-phonon coupling constant, α
- Chemical potential, μ
- Electron momentum, *k*
- Number of MC sweeps
- Number of diagram updates per sweep

## Potential Outputs
Observables:
- Total Green's function of the polaron
- Energy
- Quasi-particle weight
- Effective polaron mass

Errors:


## Packages
Feynman's variational solution to calculate polaron mobility - https://github.com/jarvist/PolaronMobility.jl

## Bibliography

1. Greitemann, J.; Pollet, L. Lecture Notes on Diagrammatic Monte Carlo for the Fröhlich Polaron. SciPost Phys. Lect. Notes 2018, 2. https://doi.org/10.21468/SciPostPhysLectNotes.2.
2. Hahn, T.; Klimin, S.; Tempere, J.; Devreese, J. T.; Franchini, C. Diagrammatic Monte Carlo Study of Fröhlich Polaron Dispersion in Two and Three Dimensions. Phys. Rev. B 2018, 97 (13), 134305. https://doi.org/10.1103/PhysRevB.97.134305.
3. Mishchenko, A. S.; Prokof’ev, N. V.; Sakamoto, A.; Svistunov, B. V. Diagrammatic Quantum Monte Carlo Study of the Fröhlich Polaron. Phys. Rev. B 2000, 62 (10), 6317–6336. https://doi.org/10.1103/PhysRevB.62.6317.
