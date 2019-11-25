# MCR

An R package for Model Reliance and Model Class Reliance (MCR, see [Fisher et al., 2019](https://arxiv.org/abs/1801.01489)). 

Given a set of models of interest (a model class), this package finds the models that rely as much, or as little as possible on variables of interest, while still performing well.

### Installation

This package uses the qp1qc package found below, to solve (possibly non-convex) quadratic programs with 1 quadratic constraint. Before using the MCR package, please run the following command.

```{r}
library(devtools)
install_github('aaronjfisher/qp1qc')
```

The MCR package can then be installed using

```{r}
library(devtools)
install_github('aaronjfisher/mcr')
```

