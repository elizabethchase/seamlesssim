#' Function to simulate a two-stage CRM
#'
#' This is a general purpose simulator for the CRM as used by the function twostage_simulator.
#' If desired, users can simulate a two-stage CRM (two rounds of toxicity assignment), rather
#' than a one-stage CRM.
#'
#' @param n_sim How many simulated trials to conduct? (positive integer)
#' @param titecrm_args This is a named list providing all of the arguments that the function
#' titesim_ss expects. See titesim_ss documentation for more information on these required components.
#' If the arguments prior, x0, and scale vary between simulations, then these arguments should be specified
#' separately using sim_specific_prior, sim_specific_x0, and sim_specific_scale, respectively.
#' @param second_stage_start_after If simulating a two-stage CRM, this positive integer indicates after how
#' many patients should we switch to the second stage. If a one-stage CRM is desired, second_stage_start_after
#' should equal n in titecrm_args.
#' @param first_stage_label A numeric value to be appended to all patients who belonged to the first stage.
#' @param sim_specific_start_id A positive integer vector containing the starting subject id for
#' each simulated trial. If provided, it must be as long as n_sim.
#' @param sim_specific_prior If provided, this is a positive numeric matrix with number of rows
#' equal to n_sim and number of columns equal to the number of dose levels, i.e. length(titecrm_args$PI).
#' Each row gives the trial-specific skeleton to use. If not provided, then a common skeleton is used,
#' taken from the value of titecrm_args$prior.
#' @param sim_specific_x0 If provided, this is a non-negative integer vector with length equal to n_sim
#' giving the starting dose level for each trial. If not provided, then a common starting dose is used,
#' taken from the value of titecrm_args$x0.
#' @param sim_specific_scale If provided, this is a positive numeric vector with length equal to
#' n_sim giving the trial-specific value of the prior scale for beta in the power model
#' p = skeleton ^ exp(beta). If not provided, then a common scale is used across simulations, taken from
#' titecrm_args$scale.
#' @param seed A positive integer seed for use prior to starting the simulations.
#' @return The function returns the a named list containing:
#' \describe{
#'   \item{all_results}{A matrix giving the individual patient outcomes from all simulated trials.}
#'   \item{first_stage_estMTD}{A vector as long as the number of simulated trials giving an integer
#' value corresponding to the estimated MTD from each trial as of the end of the first stage. A value of 0
#' indicates that all dose levels were estimated to be unsafe.}
#'   \item{second_stage_estMTD}{The same result as first_stage_estMTD but at the end of the
#' second stage.}
#'   \item{first_stage_enrollment}{A vector as long as the number of
#' simulated trials, giving the actual enrollment for that trial up to the end of the first stage.}
#'   \item{second_stage_enrollment}{A vector as long as the number of
#' simulated trials, giving the actual enrollment for that trial up to the end of the second stage. So
#' the second stage result will include the enrollment from the first stage (and equal it if that
#' simulated trial stopped in the first stage).}
#'   \item{seed}{The seed that was used by the function.}
#' }
#'
#' @references
#'
#' \insertRef{boonstra2020}{seamlesssim}
#'
#' \insertRef{dfcrm2019}{seamlesssim}
#'
#' @importFrom dplyr summarise near %>% group_by filter select
#' @importFrom data.table first last
#' @export
sim_twostage_crm = function(n_sim,
                            titecrm_args,
                            second_stage_start_after = Inf,
                            first_stage_label = 1,
                            sim_specific_start_id = NULL,
                            sim_specific_prior = NULL,
                            sim_specific_x0 = NULL,
                            sim_specific_scale = NULL,
                            seed = sample(.Machine$integer.max,1)) {

  sim_id <- NULL
  subj_id <- NULL

  beta_ests_final =
    beta_ests_end_first_stage =
    numeric(n_sim);

  index = 0;
  all_results = matrix(0, nrow = n_sim * titecrm_args$n, ncol = 6);
  colnames(all_results) = c("sim_id","subj_id","dose_num","tox_prob","tox","stage");
  # A finite-valued 'second_stage_start_after' is not well-defined if it exceeds 'titecrm_args$n'
  if(second_stage_start_after >= titecrm_args$n) {
    second_stage_start_after = Inf;
  }
  if(is.null(sim_specific_start_id)) {
    sim_specific_start_id = rep(1,n_sim);
  }
  if(is.null(sim_specific_prior)) {
    sim_specific_prior = matrix(titecrm_args$prior, nrow = n_sim, ncol = length(titecrm_args$prior), byrow = T);
  }
  if(is.null(sim_specific_x0)) {
    sim_specific_x0 = rep(titecrm_args$x0,n_sim);
  }
  if(is.null(sim_specific_scale)) {
    sim_specific_scale = rep(titecrm_args$scale,n_sim);
  }
  stopifnot(near(length(sim_specific_start_id), n_sim));
  stopifnot(near(nrow(sim_specific_prior), n_sim) &&
              near(ncol(sim_specific_prior), length(titecrm_args$PI)) &&
              (min(sim_specific_prior) >= 0) &&
              (min(sim_specific_prior) <= 1));
  stopifnot(near(length(sim_specific_x0), n_sim) &&
              all((sim_specific_x0 >= 0) &
                    (sim_specific_x0 <= length(titecrm_args$PI))));
  stopifnot(near(length(sim_specific_scale), n_sim) &&
              all(sim_specific_scale > 0));

  set.seed(seed);
  for(i in 1:n_sim) {
    #populate the trial-specific values
    curr_seed = sample(.Machine$integer.max,1);
    titecrm_args$seed = curr_seed;
    titecrm_args$prior = sim_specific_prior[i,];
    titecrm_args$x0 = sim_specific_x0[i];
    titecrm_args$scale = sim_specific_scale[i];

    #Only need to simulate trials that start at a strictly positive dose level
    if(sim_specific_x0[i] > 0) {
      big_crm = do.call(titesim_ss, args = titecrm_args)

      new_dat = cbind(big_crm$last_sim$level,
                      (big_crm$last_sim$PI)[big_crm$last_sim$level],
                      big_crm$last_sim$tox);
      colnames(new_dat) = c("level","dose","tox");
      patient_id_seq = seq_len(nrow(new_dat)) + sim_specific_start_id[i] - 1;
      index = max(index) + seq_len(nrow(new_dat));

      if(max(patient_id_seq) <= second_stage_start_after) {
        stage = first_stage_label;
      } else {
        stage = c(rep(first_stage_label, second_stage_start_after),
                  rep(first_stage_label + 1, max(patient_id_seq) - second_stage_start_after));
      }

      all_results[index,] = cbind(i,patient_id_seq,new_dat,stage);
    } else {
      #This handles the case in which the starting dose level for a simulation is 0
      big_crm = list(last_sim = list(stop.for.tox = 0,
                                     final.est = -Inf))
    }

    if(big_crm$last_sim$stop.for.tox > 0) {
      beta_ests_final[i] = -Inf;
      if(big_crm$last_sim$stop.for.tox <= second_stage_start_after) {
        #If the trial stopped for toxicity before the start of the second stage, the the estimated value of beta at the end of
        #stage 1 should be -Inf
        beta_ests_end_first_stage[i] = -Inf;
      } else {
        #Otherwise, we did not yet know that we would be stopping for toxicity during the second stage
        beta_ests_end_first_stage[i] = big_crm$last_sim$beta.hat[second_stage_start_after + 1];
      }
    } else {
      #See lines 9-11 of this function above: if second_stage_start_after is less than Inf, it's also less than titecrm_args$n
      if(second_stage_start_after < Inf) {
        beta_ests_end_first_stage[i] = big_crm$last_sim$beta.hat[second_stage_start_after + 1];
      } else {
        beta_ests_end_first_stage[i] = big_crm$last_sim$final.est;
      }
      beta_ests_final[i] = big_crm$last_sim$final.est;
    }
    rm(big_crm);
  }
  if(max(index) > 0) {
    all_results = all_results[1:max(index),,drop = F];
    enrollment_final =
      all_results %>%
      as.data.frame() %>%
      group_by(sim_id,stage) %>%
      filter(near(subj_id, last(subj_id)) | near(subj_id, first(subj_id))) %>%
      select(sim_id, stage, subj_id) %>%
      summarise(subj_id = diff(range(subj_id)) + 1) %>%
      as.data.frame();
  } else {
    all_results =
      enrollment_final =
      all_results[0,] %>%
      as.data.frame();
  }
  ptox_final = sim_specific_prior^exp(beta_ests_final);
  if(n_sim > 1) {
    #Calculate MTDs, equal to zero if stopped for toxicity
    estMTD = apply(cbind(1,abs(ptox_final - titecrm_args[["target"]]) / (ptox_final <= (titecrm_args[["target"]] + titecrm_args[["no.exceed"]]))),1,which.min)-1;
    estMTD[beta_ests_end_first_stage == -Inf] = 0;
  } else {
    if(beta_ests_end_first_stage == -Inf) {
      estMTD = 0;
    } else {
      estMTD = which.min(c(1,abs(ptox_final - titecrm_args[["target"]]) / (ptox_final <= (titecrm_args[["target"]] + titecrm_args[["no.exceed"]])))) - 1;
    }
  }

  first_stage_enrollment = numeric(n_sim);
  first_stage_enrollment[filter(enrollment_final,near(stage, first_stage_label))[,"sim_id"]] =
    filter(enrollment_final,near(stage, first_stage_label))[,"subj_id"]
  second_stage_enrollment = rep(0, n_sim);

  if(second_stage_start_after == Inf) {

    first_stage_estMTD = estMTD;
    second_stage_estMTD = rep(NA, n_sim);

  } else {
    if(any(beta_ests_end_first_stage > (-Inf))) {
      second_stage_enrollment[which(beta_ests_end_first_stage > (-Inf))] =
        filter(enrollment_final,near(stage, 2))[,"subj_id"];
    }

    second_stage_estMTD = estMTD;
    ptox_end_first_stage = sim_specific_prior^exp(beta_ests_end_first_stage);
    if(n_sim > 1) {
      #Calculate MTDs, equal to zero if stopped for toxicity
      estMTD = apply(cbind(1,abs(ptox_end_first_stage - titecrm_args[["target"]]) / (ptox_end_first_stage <= (titecrm_args[["target"]] + titecrm_args[["no.exceed"]]))),1,which.min)-1;
      estMTD[beta_ests_end_first_stage == -Inf] = 0;
    } else {
      if(beta_ests_end_first_stage == -Inf) {
        estMTD = 0;
      } else {
        estMTD = which.min(c(1,abs(ptox_end_first_stage - titecrm_args[["target"]]) / (ptox_end_first_stage <= (titecrm_args[["target"]] + titecrm_args[["no.exceed"]])))) - 1;
      }
    }
    first_stage_estMTD = estMTD;
  }

  list(all_results = all_results,
       first_stage_estMTD = first_stage_estMTD,
       second_stage_estMTD = second_stage_estMTD,
       first_stage_enrollment = first_stage_enrollment,
       second_stage_enrollment = second_stage_enrollment,
       seed = seed);
}


count_stan_divergences = function(stan_fit, sum_chains = TRUE) {
  foo = get_sampler_params(stan_fit, inc_warmup = FALSE);
  n_draws = lapply(foo, nrow)[[1]];
  sum(unlist(lapply(foo,"[",i = 1:n_draws, j = "divergent__")));
}
