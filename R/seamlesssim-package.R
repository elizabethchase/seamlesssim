#' The 'seamlesssim' package.
#'
#' @description seamlesssim simulates complex seamless Phase I/II oncology trials
#' as discussed in the article by Boonstra et al. (2021). It allows clinical
#' trialists to determine operating characteristics of trials that assess both
#' toxicity and efficacy with a range of different design and analytic approaches.
#' For more information, see Boonstra et al. and the vignette.
#'
#' @docType package
#' @name seamlesssim-package
#' @aliases seamlesssim
#' @useDynLib seamlesssim, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#'
#' \insertRef{boonstra2020}{seamlesssim}
#'
#' \insertRef{boonstra2020b}{seamlesssim}
#'
#' \insertRef{dfcrm2019}{seamlesssim}
#'
#' \insertRef{rstan2020}{seamlesssim}
NULL
