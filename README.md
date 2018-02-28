# The energetics of fish growth and how it constrains food-web trophic structure

This repository contains code and data needed to reproduce the article:

**Barneche DR, Allen AP** (in press) The energetics of fish growth and how it constrains food-web trophic structure. *Ecology Letters*. doi: 10.1111/ele.12947.  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1186364.svg)](https://doi.org/10.5281/zenodo.1186364)  

## Instructions

All analyses were done in `R`. To compile the paper, including figures and tables we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies=TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies=TRUE)
```

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

And then install the package `fontcm`, via `extrafont`. This installs the font `CM Roman` we use in our figures (for more information on see [these instructions](https://cran.r-project.org/web/packages/fontcm/README.html):

```r
extrafont::font_install('fontcm')
```

Then, to generate all figures, analyses, and manuscript (.docx, using Rmarkdown), simply do:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that the Bayesian analysis in this paper might take a couple of days to run on a regular computer.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename="build.R")
```

### This paper was produced using the following software and associated packages:
```
R version 3.4.3 (2017-11-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-2 extrafont_0.17     maps_3.2.0         fontcm_1.1         LoLinR_0.0.0.9000  loo_1.1.0         
 [7] brms_2.1.0         Rcpp_0.12.15       rstan_2.17.3       StanHeaders_2.17.2 ggplot2_2.2.1      plyr_1.8.4        

loaded via a namespace (and not attached):
 [1] mvtnorm_1.0-7        lattice_0.20-35      gtools_3.5.0         zoo_1.8-1            lmtest_0.9-35       
 [6] assertthat_0.2.0     digest_0.6.15        mime_0.5             R6_2.2.2             stats4_3.4.3        
[11] coda_0.19-1          colourpicker_1.0     pillar_1.1.0         rlang_0.1.6          lazyeval_0.2.1      
[16] miniUI_0.1.1         extrafontdb_1.0      Matrix_1.2-12        DT_0.4               shinythemes_1.1.1   
[21] shinyjs_1.0          stringr_1.2.0        htmlwidgets_1.0      igraph_1.1.2         munsell_0.4.3       
[26] shiny_1.0.5          compiler_3.4.3       httpuv_1.3.5         pkgconfig_2.0.1      base64enc_0.1-3     
[31] rstantools_1.4.0     htmltools_0.3.6      tibble_1.4.2         gridExtra_2.3        threejs_0.3.1       
[36] matrixStats_0.53.0   remake_0.3.0         crayon_1.3.4         dplyr_0.7.4          grid_3.4.3          
[41] nlme_3.1-131         xtable_1.8-2         Rttf2pt1_1.3.5       gtable_0.2.0         magrittr_1.5        
[46] storr_1.1.3          scales_0.5.0         stringi_1.1.6        reshape2_1.4.3       bindrcpp_0.2        
[51] dygraphs_1.1.1.4     xts_0.10-1           tools_3.4.3          glue_1.2.0           markdown_0.8        
[56] shinystan_2.4.0      crosstalk_1.0.0      rsconnect_0.8.5      abind_1.4-5          parallel_3.4.3      
[61] yaml_2.1.16          inline_0.3.14        colorspace_1.3-2     bridgesampling_0.4-0 bayesplot_1.4.0     
[66] bindr_0.1            Brobdingnag_1.2-4   
```

## How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  

## Bug reporting
* Please [report any issues or bugs](https://github.com/dbarneche/FishGrowth/issues).
