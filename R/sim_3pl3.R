#' Function to simulate a 3+3 trial design
#'
#' This function is an efficient simulator of the 3+3 design. It is intended to be called from
#' twostage_simulator rather than by the user directly.
#'
#' @param n_sim How many simulated trials to conduct? (positive integer)
#' @param true_tox_curve A positive numeric vector that contains the true generating dose-toxicity
#' curve to simulate the data from. This also gives the number of dose levels.
#' @param stage_label A numeric value that can be arbitrary but is intended to take on integer
#' values equal to either 1 or 2, corresponding to the stage of the trial.
#' @param sim_specific_start_id A positive integer vector that contains the starting subject id
#' for each simulated trial. If provided, it must be as long as 'n_sim.' If left blank, the default
#' starting patient is patient 1.
#' @param sim_specific_dose_start A positive integer vector that contains the starting dose level
#' for each simulated trial. If provided, it must be as long as 'n_sim' and take on values between
#' 0 (indicating that no patients should be enrolled) to length(true_tox_curve). If left blank, the
#' default starting dose is dose 1.
#' @param seed A positive integer seed for use prior to starting the simulations.
#' @return A named list with entries:
#' \describe{
#'    \item{all_results}{A matrix giving the individual patient outcomes from all simulated trials.}
#'    \item{estMTD}{A vector as long as the number of simulated trials giving an integer value
#' corresponding to the estimated MTD from each trial; a value of 0 indicates that all dose
#' levels were estimated to be unsafe.}
#'    \item{enrollment}{A vector as long as the number of simulated trials giving the total enrollment
#' for that trial.}
#'    \item{seed}{The seed that was used by the function.}
#' }
#'
#' @references
#' \insertRef{boonstra2020}{seamlesssim}
#'
#'
#' @importFrom dplyr near
#' @importFrom stats binomial
#' @export
sim_3pl3 = function(n_sim,
                    true_tox_curve,
                    stage_label = 1,
                    sim_specific_start_id = NULL,
                    sim_specific_dose_start = NULL,
                    seed = sample(.Machine$integer.max,1)) {

  n_dose = length(true_tox_curve);
  #Check that arguments are equal length
  if(is.null(sim_specific_start_id)) {
    sim_specific_start_id = rep(1, n_sim);
  }
  if(is.null(sim_specific_dose_start)) {
    sim_specific_dose_start = rep(1, n_sim);
  }
  stopifnot(near(length(sim_specific_dose_start), n_sim));
  stopifnot(near(length(sim_specific_start_id), n_sim));
  stopifnot(min(sim_specific_dose_start) >= 0 && max(sim_specific_dose_start) <= n_dose);

  #Update n_sim and note which sim_ids actually have to be simulated (any that start at dose level '0' need not be)
  estMTD = sim_specific_dose_start;
  enrollment = numeric(n_sim);
  sim_id = which(sim_specific_dose_start > 0);
  sim_specific_start_id = sim_specific_start_id[sim_id];
  sim_specific_dose_start = sim_specific_dose_start[sim_id];
  n_sim = length(sim_id);

  set.seed(seed);
  #Now create all possible dose assignments a priori
  all_results = matrix(0, nrow = 6 * n_dose * n_sim,ncol = 7);
  colnames(all_results) = c("sim_id","subj_id","dose_num","tox_prob","tox","stage","starting_dose");

  all_results[,"sim_id"] = rep(sim_id,each = 6 * n_dose);
  all_results[,"dose_num"] = rep(rep(1:n_dose, times = n_sim), each = 6);
  all_results[,"tox_prob"] = true_tox_curve[all_results[,"dose_num"]];
  all_results[,"tox"] = rbinom(nrow(all_results),1,all_results[,"tox_prob"]);
  all_results[,"stage"] = stage_label;
  all_results[,"starting_dose"] =
    near(all_results[,"dose_num"], rep(sim_specific_dose_start, each = 6 * n_dose));
  all_results = all_results[order(all_results[,"sim_id"],-all_results[,"starting_dose"],all_results[,"subj_id"]),1:6];
  reorder_all_rows = rep(NA,nrow(all_results));

  for(i in 1:n_sim) {
    sim_index = which(near(all_results[,"sim_id"], sim_id[i]));
    curr_results = all_results[sim_index,];
    curr_dose = sim_specific_dose_start[i];
    curr_dose_index = row_order = 1:3;
    max_safe_dose = n_dose;
    #Keep track of whether six have already been at this dose or not
    six_at_dose = rep(F,n_dose);
    finished = F;
    while(!finished) {
      #If no toxicity at all
      if(sum(curr_results[curr_dose_index,"tox"]) == 0) {
        #escalate if the next dose is safe
        if(curr_dose < max_safe_dose) {
          curr_dose = curr_dose + 1;
          curr_dose_index = which(curr_results[,"dose_num"] == curr_dose)[1:3];
          row_order = c(row_order,curr_dose_index);
          next;
        } else {
          #if six are at the current dose, this is the MTD
          if(six_at_dose[curr_dose]) {
            finished = T;
            next;
            #otherwise enroll three more
          } else {
            curr_dose_index = which(curr_results[,"dose_num"] == curr_dose);
            row_order = c(row_order,curr_dose_index[4:6]);
            six_at_dose[curr_dose]=T;
            next;
          }
        }
        #If one toxicity
      } else if(sum(curr_results[curr_dose_index,"tox"]) == 1) {
        #If already six at the current dose
        if(six_at_dose[curr_dose]) {
          #escalate if safe to
          if(curr_dose<max_safe_dose) {
            curr_dose = curr_dose + 1;
            curr_dose_index = which(curr_results[,"dose_num"] == curr_dose)[1:3];
            row_order = c(row_order,curr_dose_index);
            next;
            #otherwise this is the MTD
          } else {
            finished = T;
            next;
          }
          #Otherwise enroll three more
        } else {
          curr_dose_index = which(curr_results[,"dose_num"] == curr_dose);
          row_order = c(row_order,curr_dose_index[4:6]);
          six_at_dose[curr_dose]=T;
          next;
        }
        #Otherwise there are 2 or more toxicities, so de-escalate
      } else {
        curr_dose = max_safe_dose = curr_dose - 1;
        #stop if everything is too toxic
        if(curr_dose == 0) {
          finished = T;
          next;
        } else {
          #if six are already at the lower dose, this is the MTD
          if(six_at_dose[curr_dose]) {
            finished = T;
            next;
          } else {
            #if you have already visited this dose with three patients
            if(curr_dose%in%curr_results[row_order,"dose_num"]) {
              curr_dose_index = which(curr_results[,"dose_num"] == curr_dose);
              row_order = c(row_order,curr_dose_index[4:6]);
              six_at_dose[curr_dose]=T;
              next;
            } else {
              #otherwise you havent ever visited this dose
              curr_dose_index = which(curr_results[,"dose_num"] == curr_dose)[1:3];
              row_order = c(row_order,curr_dose_index);
              next;
            }
          }
        }
      }
    }
    estMTD[sim_id[i]] = curr_dose;
    reorder_all_rows[min(sim_index)-1+(1:length(row_order))] = min(sim_index)-1 + row_order;
  }
  all_results = all_results[reorder_all_rows[which(!is.na(reorder_all_rows))],];

  enrollment = as.numeric(tapply(all_results[,"sim_id"],all_results[,"sim_id"],length))

  all_results[,"subj_id"] = unlist(mapply(rep, x = sim_specific_start_id, times = enrollment)) - 1 +
    unlist(sapply(diff(c(0,which(!duplicated(all_results[,"sim_id"],fromLast=T)))),seq,from=1,simplify=F))

  return(list(all_results = all_results,
              estMTD = estMTD,
              enrollment = enrollment,
              seed = seed));
}
