# mclustAddonsEHA: A fork of CRAN mclustAddons with modifications required at EcoHealth Alliance

- Renamed to avoid compatibility issues
- Changes will be offered upstream to mclust and mclustAddons
- Install from 

    install.packages('mclustAddonsEHA', repos = c('https://ecohealthalliance.r-universe.dev', getOptions('repos')))
 
- Key features:
    - Incorporation of priors into densityMclustBounded
---

# mclustAddons

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/mclustAddons)](https://cran.r-project.org/package=mclustAddons)
[![CRAN\_MonthlyDownloads](http://cranlogs.r-pkg.org/badges/mclustAddons)](https://cran.r-project.org/package=mclustAddons)

An [R](https://www.r-project.org/) package extending the functionality of the [mclust](https://mclust-org.github.io/mclust/index.html) package (Scrucca et al. 2016) by including:

- density estimation for data with bounded support using a transform-based approach to Gaussian mixture density estimation (Scrucca, 2019);

- modal clustering using modal EM algorithm for Gaussian mixtures (Scrucca, 2021);

- entropy estimation via Gaussian mixture modeling (Robin & Scrucca, 2023).
  
## Installation

You can install the released version of **mclustAddons** from CRAN using:

```{r}
install.packages("mclustAddons")
```

## Usage

For an introduction to the main functions and several examples see the vignette **A quick tour of mclustAddons**, which is available as

```{r}
vignette("mclustAddons")
```

The vignette is also available in the *Vignette* section on the navigation bar on top of the package's web page.

<br><br>

## References

Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, *The R Journal*, 8/1, 205-233. https://doi.org/10.32614/RJ-2016-021

Scrucca L. (2019) A transformation-based approach to Gaussian mixture density estimation for bounded data, *Biometrical Journal*, 61:4, 873–888. https://doi.org/10.1002/bimj.201800174

Scrucca L. (2021) A fast and efficient Modal EM algorithm for Gaussian mixtures. *Statistical Analysis and Data Mining*, 14:4, 305–314. https://doi.org/10.1002/sam.11527

Robin S. and Scrucca L. (2023) Mixture-based estimation of entropy. *Computational Statistics & Data Analysis*, 177, 107582. https://doi.org/10.1016/j.csda.2022.107582
