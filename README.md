Delta-statistics
================
Kalle Leppälä

This repository contains the code for the article (add a link here).

The article is joint work between Kalle Leppälä, Flavio Augusto da Silva
Coelho, Michaela Richter, Victor A. Albert and Charlotte Lindqvist.

## Simulation

- [tools.md](tools.md) has the R functions and documentation used in all
  the coalescence simulations.
- [S_big.md](S_big.md) is the simulation and documentation in tree $S$
  using 1,000,000 allelic patterns.
- [S_medium.md](S_medium.md) is the simulation in tree $S$ using 100,000
  allelic patterns.
- [S_small.md](S_small.md) is the simulation in tree $S$ using 10,000
  allelic patterns.
- [A_big.md](A_big.md) is the simulation and documentation in tree $A$
  using 1,000,000 allelic patterns.
- [A_medium.md](A_medium.md) is the simulation in tree $A$ using 100,000
  allelic patterns.
- [A_small.md](A_small.md) is the simulation in tree $A$ using 10,000
  allelic patterns.
- [Q_big.md](S_big.md) is the simulation and documentation in tree $Q$
  using 1,000,000 allelic patterns.
- [Q_medium.md](S_medium.md) is the simulation in tree $Q$ using 100,000
  allelic patterns.
- [Q_small.md](S_small.md) is the simulation in tree $Q$ using 10,000
  allelic patterns.
- [DFOIL.md](DFOIL.md) is the simulation in tree $S$ using
  $D_\text{FOIL}$.
- [dependence.md](dependence.md) explores the dependence of the
  $D_\text{FOIL}$- and $\Delta$-statistics on the synchronization
  assumption (2).

## Empirical example

- [Delta.py](Delta.py) contains Python code for computing
  $\Delta$-statistics form PLINK `.traw` -files, either in windows or
  using the block jackknife.
- [Summary_script_without_singletons.R](Summary_script_without_singletons.R)
  contains R code used for classifying the $\Delta$-statistic
  signatures, without singleton patterns.
- [Summary_script_singletons.R](Summary_script_singletons.R) contains R
  code used for classifying the $\Delta$-statistic signatures, with
  singleton patterns.
