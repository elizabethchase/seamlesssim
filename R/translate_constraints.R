#' Function to translate prior isotonic regression constraints into isotonic regression model parameters
#'
#' This function translates a desired value for the prior mean and effective prior sample size of a
#' set of probabilities that are assumed to be monotonically increasing in some arbitrary way into 
#' a vector of shape and rate parameters such that each individual probability is marginally 
#' distributed as a beta random variable and the probabilities are ordered in an increasing fashion.
#' 
#' @param prior_mean A numeric vector comprised of increasing values between 0 and 1. This is a vector 
#' giving the desired prior means of the efficacy probabilities. 
#' @param prior_n_per A positive numeric scalar giving the effective prior sample size. 

#' @return The function returns vectors of the shape and rate parameters of gamma distributions that the 
#' function "bayesian_isotonic()" requires. Literally, what is returned is a named list:
#' \describe{
#'   \item{gamma_shape}{The set of shape parameters to provide to the underlying gamma distribution}
#'   \item{gamma_rate}{The set of rate parameters. In order for the set of means (mu_j) to be both 
#'   marginally distributed as beta random variables and isotonically ordered, gamma_rate must be
#'   identically equal to 1.}
#' }
#' @export
translate_constraints = function(prior_mean, prior_n_per = 1) {
  if(any(diff(prior_mean) <= 0) ||
     min(prior_mean) < 0 || 
     max(prior_mean) > 1 ||
     length(prior_n_per) > 1 ||
     prior_n_per <= 0) {
    stop("'prior_mean' must be an increasing vector in [0,1] and 'prior_n_per' must be a positive scalar");
  }
  
  list(gamma_shape = diff(c(0,prior_mean,1)) * prior_n_per, 
       gamma_rate = 1 + numeric(length(prior_mean) + 1));
}