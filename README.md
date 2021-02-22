# MultiBaC R package repository

[![](https://img.shields.io/badge/bioc%20release-6.12.2-green.svg)](https://www.bioconductor.org/packages/MultiBaC)
[![download](http://www.bioconductor.org/shields/downloads/release/MultiBaC.svg)](https://bioconductor.org/packages/stats/bioc/MultiBaC)
[![](https://img.shields.io/github/last-commit/ConesaLab/MultiBaC.svg)](https://github.com/ConesaLab/MultiBaC/commits/master)
[![dependencies](http://bioconductor.org/shields/dependencies/release/MultiBaC.svg)](http://bioconductor.org/packages/release/bioc/html/MultiBaC.html#since)

This repository contains the `R` package [now hosted on
Bioconductor](http://bioconductor.org/packages/release/bioc/html/MultiBaC.html)
and our current `GitHub` version. A comprehensive user's guide can be found at [http://www.bioconductor.org/packages/release/bioc/vignettes/MultiBaC/inst/doc/MultiBaC.html](http://www.bioconductor.org/packages/release/bioc/vignettes/MultiBaC/inst/doc/MultiBaC.html)

## Installation

**(Mac OS Users Only:)** Ensure you have installed
[XQuartz](https://www.xquartz.org/) first.

Make sure you have the latest R version and the latest `BiocManager`
package installed following [these
instructions](https://www.bioconductor.org/install/) (if you use legacy
R versions (\<=3.5.0) refer to the instructions at the end of the
mentioned page).

``` r
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
## ensure the following returns TRUE, or follow guidelines
BiocManager::valid()
```

#### Latest Bioconductor Release

You can then install `MultiBaC` using the following
code:

```r
## install mixOmics
BiocManager::install('mixOmics')
```

#### `GitHub` Versions

MultiBaC can be also directly installed from this repository:

```r
install.packages("devtools")
evtools::install_bitbucket("ConesaLab/MultiBaC")
```

## Contribution

#### Bug reports and pull requests

To report a bug (or offer a solution for a bug\!):
<https://github.com/ConesaLab/MultiBaC/issues>. We fully welcome and
appreciate well-formatted and detailed pull requests. Preferrably with
tests on our datasets.

## Citation

[1] Ugidos, M., Tarazona, S., Prats-Montalbán, J. M., Ferrer, A., & Conesa, A. (2020). MultiBaC: A strategy to remove batch effects between different omic data types. Statistical Methods in Medical Research. https://doi.org/10.1177/0962280220907365
[2] Nueda MJ, Ferrer A, Conesa A. ARSyN: A method for the identification and removal of systematic noise in multifactorial time course microarray experiments. Biostatistics. 2012;13:553–66.