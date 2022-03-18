# Hitman

[![Build Status](https://travis-ci.com/jdreyf/Hitman.svg?branch=master)](https://travis-ci.com/jdreyf/Hitman)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/Hitman/master.svg)](https://codecov.io/github/jdreyf/Hitman?branch=master)

## Install
On Windows, you should have [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Install `Hitman` from GitHub using `remotes` within R. You must install `remotes`, e.g. with `install.packages("remotes")`, if you haven't before. `Hitman` depends on `ezlimma` which depends on `limma` so you must also install these using instruction below if you haven't before.
```
#if haven't already installed limma
install.packages("BiocManager") #if haven't already installed BiocManager
library(BiocManager)
BiocManager::install("limma")

library(remotes)
remotes::install_github(repo="jdreyf/ezlimma", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
remotes::install_github(repo="jdreyf/Hitman", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

## Usage
The vignette presents a tutorial. To see the vignette:
```
library(limma)
library(ezlimma)
library(Hitman)
library(rmarkdown)
browseVignettes(package="Hitman")
```
and click on HTML.
