#' The 'parTimeROC' package.
#'
#' @description The goal of parTimeROC is to store methods and procedures needed to run the time-dependent ROC analysis parametrically. This package adopts two different theoretical framework to produce the ROC curve which are from the proportional hazard model and copula function. Currently, this package only able to run analysis for single covariate/biomarker with survival time. The future direction for this work is to be able to include analysis for multiple biomarkers with longitudinal measurements.
#'
#' @name parTimeROC-package
#' @aliases parTimeROC
#' @useDynLib parTimeROC, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.3. https://mc-stan.org
#'
"_PACKAGE"
