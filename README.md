# Hitman

[![Build Status](https://travis-ci.org/jdreyf/Hitman.svg?branch=master)](https://travis-ci.org/jdreyf/Hitman)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/Hitman/master.svg)](https://codecov.io/github/jdreyf/Hitman?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## Install
On Windows, you should have [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Install `Hitman` from GitHub using `remotes` within R. You must install `remotes`, e.g. with `install.packages("remotes")`, if you haven't before. `Hitman` depends on `limma` so you must also install this using instruction below if you haven't before.
```
#if haven't already installed limma
install.packages("BiocManager") #if haven't already installed BiocManager
library(BiocManager)
BiocManager::install("limma")

library(remotes)
remotes::install_github(repo="jdreyf/Hitman", build_opts = c("--no-resave-data", "--no-manual"))
```

## Usage
The vignette presents a tutorial. To see the vignette:
```
library(limma)
library(Hitman)
library(rmarkdown)
browseVignettes(package="Hitman")
```
and click on HTML.

## Code of Conduct

This project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.