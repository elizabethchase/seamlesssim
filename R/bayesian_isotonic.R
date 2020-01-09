#' Function to fit a Bayesian isotonic regression
#'
#' This function fits a Bayesian isotonic regression: it calculates the posterior distribution
#' of a set of probabilities that are order constrained but each a priori marginally
#' distributed as beta random variables.
#'
#' @param data_grouped A data.frame or tibble that contain columns x (used as category
#' labels), y (the number of successes or events), and n (the number of trials)
#' @param stan_args A named list of arguments to pass to the stan function. It should include the
#' named elements:
#' \describe{
#'     \item{local_dof_stan}{This has default value 1.}
#'     \item{global_dof_stan}{This has default value 1.}
#'     \item{alpha_scale_stan}{This has default value 1.}
#'     \item{slab_precision_stan}{This has default value 1.}
#' }
#' For more information on these arguments, please see rstan documentation.
#' @param sample_from_prior_only A logical value. If TRUE, then the provided values of 'x'
#' and 'n' will be ignored and draws will only be sampled from the prior
#' distribution. Defaults to FALSE.
#' @param conf_level Numeric, (0, 1). This level credible interval will be
#' returned (in data_grouped) based upon the empirical quantiles. The default is 0.5.
#' @param conf_level_direction This is a string equal to 'both', 'upper',
#' or 'lower' indicating the desired direction of the credible interval. The default is 'both'.
#' @param verbose A logical value. If TRUE, then all posterior draws will be returned
#' by the function. Defaults to FALSE.
#' @param n_mc_warmup See documentation for twostage_simulator. The default here is 1000.
#' @param n_mc_samps See documentation for twostage_simulator. The default here is 2000.
#' @param mc_chains See documentation for twostage_simulator. The default here is 4.
#' @param mc_thin See documentation for twostage_simulator. The default here is 1.
#' @param mc_stepsize See documentation for twostage_simulator. The default here is 0.1.
#' @param mc_adapt_delta See documentation for twostage_simulator. The default here is 0.8.
#' @param mc_max_treedepth See documentation for twostage_simulator. The default here is 15.
#' @param ntries See documentation for twostage_simulator. The default here is 2.
#' @param return_as_stan_object A logical value. If TRUE, then the function returns an
#' object of class 'stanfit'. If FALSE, then a summary of results will be
#' returned. Defaults to FALSE.
#' @param tol A small positive number (double). The default is the square root of the machine precision.
#' @return If return_as_stan_object = TRUE, then an object of class stanfit is returned.
#' This is useful for the initial compilation of the stan model. Otherwise, the function returns the following
#' named list containing the arguments:
#' \describe{
#'    \item{data_grouped}{The provided argument of the same name but with more columns added
#' that provide various summaries of the probabilities at each group level}
#'    \item{conf_level, conf_level_direction}{Arguments provided by the user for the confidence
#'    interval}
#'    \item{stan_args}{The list of arguments that were passed to the stan function}
#'    \item{accepted_divergences}{The number of divergent transitions from the model fit the results
#' of which were actually returned}
#'    \item{max_divergences}{The maximum observed number of divergent transitions from the ntries
#' number of model fits that were attempted}
#'    \item{rhat}{The largest value of the Gelman-Rubin diagnostic across all parameters}
#'    \item{number_nan}{The number of draws that were 0/0, e.g. due to underflow}
#'    \item{all_draws}{NA if verbose == FALSE, otherwise all_draws is a matrix of draws from the
#' posterior distribution, which may be large}
#'    \item{chain_run_times_secs, total_run_time_secs}{Length of time (seconds) for chain runs and total runs}
#' }
#' @references PS Boonstra, DR Owen, and J Kang, "The isotonic horseshoe prior for modeling binary outcomes." Arxiv, 2020.
#' @importFrom dplyr pull near %>% mutate bind_cols
#' @importFrom tibble as_tibble
#' @importFrom stats quantile
#' @import rstan
#' @export
bayesian_isotonic = function(data_grouped = NULL,
                             stan_args = list(
                               local_dof_stan = 1,
                               global_dof_stan = 1,
                               alpha_scale_stan = 1,
                               slab_precision_stan = 1),
                             sample_from_prior_only = F,
                             conf_level = 0.50,
                             conf_level_direction = "both",
                             verbose = F,
                             n_mc_warmup = 1e3,
                             n_mc_samps = 2e3,
                             mc_chains = 4,
                             mc_thin = 1,
                             mc_stepsize = 0.1,
                             mc_adapt_delta = 0.8,
                             mc_max_treedepth = 15,
                             ntries = 2,
                             return_as_stan_object = F,
                             tol = .Machine$double.eps^0.5) {

  stopifnot(c("y","n") %in% colnames(data_grouped));
  y <- NULL
  stopifnot(all(pull(data_grouped,y) >= 0) &&
              all(pull(data_grouped,y) <= pull(data_grouped,n)) &&
              all(pull(data_grouped,n) >= 0) );

  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;

  while(curr_try <= ntries) {

    curr_fit <- sampling(object = stanmodels$iso_horseshoe,
                     data = c(list(n_groups_stan = nrow(data_grouped),
                                   n_per_group_stan = as.array(pull(data_grouped,n)),
                                   y_stan = as.array(pull(data_grouped,y)),
                                   only_prior_stan = as.integer(sample_from_prior_only)),
                              stan_args),
                     warmup = n_mc_warmup,
                     iter = n_mc_samps + n_mc_warmup,
                     chains = mc_chains,
                     thin = mc_thin,
                     verbose = F,
                     control = list(stepsize = mc_stepsize,
                                    adapt_delta = mc_adapt_delta,
                                    max_treedepth = mc_max_treedepth),
                     refresh=0);

    curr_divergences = count_stan_divergences(curr_fit);#unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    max_divergences = max(max_divergences,curr_divergences,na.rm=T);

    rhat_check = max(summary(curr_fit)$summary[,"Rhat"])
    # Originally, the break conditions were based upon having both no divergent
    # transitions as well as a max Rhat (i.e. gelman-rubin diagnostic)
    # sufficiently close to 1. I subsequently changed the conditions to be based
    # only upon the first, which is reflecte by setting rhat = T immediately below.
    break_conditions = c(divergence = F, rhat = T);

    if(near(curr_divergences, 0)) {#corresponds to zero divergent transitions
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      if(!return_as_stan_object) {
        chain_run_times_secs = rowSums(rstan::get_elapsed_time(curr_fit));
        total_run_time_secs = max(chain_run_times_secs);
        foo = rstan::extract(curr_fit);

        number_nan = sum(rowSums(is.na(foo$xi)) > 0);
        mean_prob = colMeans(foo$xi, na.rm = T);

        if(conf_level_direction == "both") {
          quantile_probs <- apply(foo$xi,2,quantile,p=1/2+c(-conf_level,0,conf_level)/2,na.rm = T);
        } else if(conf_level_direction == "lower") {
          quantile_probs <- rbind(apply(foo$xi,2,quantile,p=c(1-conf_level,1/2),na.rm = T),
                                  "100%" = 1);
        } else {
          quantile_probs <- rbind("0%" = 0,
                                  apply(foo$xi,2,quantile,p=c(1/2,conf_level),na.rm = T));
        }
        quantile_probs = t(quantile_probs);
        colnames(quantile_probs) = c("model_lower_ci_prob", "model_median_prob", "model_upper_ci_prob");
      } else {
        foo = rstan::extract(curr_fit);
        number_nan = sum(rowSums(is.na(foo$xi)) > 0);
      }
    }
    if(all(break_conditions)) {
      break;
    }
  }

  if(accepted_divergences > 0) {
    warning(paste0("there were ", accepted_divergences, " divergent transitions"));
  }
  if(number_nan > 0) {
    warning(paste0("there were ", number_nan, " draws in which one or more elements of xi were NaN"));
  }

  if(return_as_stan_object) {
    curr_fit;
  } else {

    data_grouped =
      data_grouped %>%
      mutate(emp_mean_prob = y/n) %>%
      bind_cols(model_mean_prob = mean_prob,
                as_tibble(quantile_probs));

    draws_delta = t(apply(cbind(0,0,foo$xi),1,diff))[,-1,drop = F];

    if(verbose) {
      c(list(data_grouped = data_grouped,
             conf_level = conf_level),
        stan_args,
        list(sample_from_prior_only = sample_from_prior_only,
             accepted_divergences = accepted_divergences,
             max_divergences = max_divergences,
             rhat = rhat_check,
             number_nan = number_nan,
             all_draws = foo,
             chain_run_times_secs = chain_run_times_secs,
             total_run_time_secs = total_run_time_secs));
    } else {
      c(list(data_grouped = data_grouped,
             conf_level = conf_level),
        stan_args,
        list(sample_from_prior_only = sample_from_prior_only,
             accepted_divergences = accepted_divergences,
             max_divergences = max_divergences,
             rhat = rhat_check,
             number_nan = number_nan,
             all_draws = NA,
             chain_run_times_secs = chain_run_times_secs,
             total_run_time_secs = total_run_time_secs));

    }
  }
}
