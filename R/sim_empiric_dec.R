#' A function to simulate a dose expansion cohort trial design
#'
#' A general purpose simulator for dose expansion cohorts (DECs). It is intended to be called
#' from twostage_simulator rather than by the user directly.
#'
#' @param n_sim How many simulated trials to conduct? (positive integer)
#' @param true_tox_curve A positive numeric vector that contains the true generating dose-toxicity
#' curve to simulate the data from. This also gives the number of dose levels.
#' @param stage_label A numeric value that can be arbitrary but is intended to take on integer
#' values equal to either 1 or 2, corresponding to the stage of the trial.
#' @param sim_specific_start_id A positive integer vector that contains the starting subject id
#' for each simulated trial. If provided, it must be as long as n_sim.
#' @param sim_specific_dose_start A positive integer vector that contains the starting dose level
#' for each simulated trial. If provided, it must be as long as n_sim and take on values between
#' 0 (indicating that no patients should be enrolled) to length(true_tox_curve).
#' @param max_n_per_dec A positive integer giving the maximum enrollment for each DEC, if no
#' stopping for toxicity occurs.
#' @param module_rule Currently the only valid choice for this is the string "local", corresponding
#' to a local stopping rule: if the empiric proportion of DLTs at the current dose level ever
#' exceeds thresh_decrease, then de-escalate.
#' @param thresh_decrease A numeric value between 0 and 1. This is the maximum tolerance for the
#' observed proportion of DLTs at any given dose level. If this threshold is ever exceeded, the
#' DEC will de-escalate if possible or stop the trial entirely if not.
#' @param first_patient_look An integer greater than or equal to 0 and less than or equal to max_n_per_dec.
#' At a given dose level, at what point should the module_rule start taking action? This is an
#' ad-hoc way to prevent a scenario such as "The first patient experienced a DLT, therefore because
#' all patients at this dose level have experienced a DLT, we should de-escalate."
#' @param seed A positive integer seed for use prior to starting the simulations.
#' @return A named list with entries:
#' \describe{
#'    \item{all_results}{A matrix giving the individual patient outcomes from all simulated trials.}
#'    \item{estMTD}{A vector as long as the number of simulated trials giving an integer value corresponding
#'    to the estimated MTD from each trial; a value of 0 indicates that all dose levels were estimated to be unsafe.}
#'    \item{enrollment}{A vector as long as the number of simulated trials giving the total enrollment for that trial.}
#'    \item{seed}{The seed that was used by the function.}
#' }
#' @importFrom stats binomial
#' @export
sim_empiric_dec = function(n_sim,
                           true_tox_curve,
                           stage_label = 1,
                           sim_specific_start_id = NULL,
                           sim_specific_dose_start = NULL,
                           max_n_per_dec,
                           module_rule = "local",
                           thresh_decrease = 1/3,
                           first_patient_look = 0,
                           seed = sample(.Machine$integer.max,1)) {

  set.seed(seed);
  n_dose = length(true_tox_curve);
  #Check that arguments are equal length
  if(is.null(sim_specific_start_id)) {
    sim_specific_start_id = rep(1, n_sim);
  }
  if(is.null(sim_specific_dose_start)) {
    sim_specific_dose_start = rep(n_dose, n_sim);
  }
  stopifnot(length(sim_specific_dose_start) == n_sim);
  stopifnot(length(sim_specific_start_id) == n_sim);
  stopifnot(min(sim_specific_dose_start) >= 0 && max(sim_specific_dose_start) <= n_dose);
  stopifnot(first_patient_look >= 0 && first_patient_look < max_n_per_dec);

  #Update n_sim and note which sim_ids actually have to be simulated (any that start at dose level '0' need not be)
  estMTD = sim_specific_dose_start;
  enrollment = numeric(n_sim);
  sim_id = which(sim_specific_dose_start > 0);
  sim_specific_start_id = sim_specific_start_id[sim_id];
  sim_specific_dose_start = sim_specific_dose_start[sim_id];
  n_sim = length(sim_id);

  #Now create all remaining possible outcomes, assuming that no de-escalation rules are in place. Pruning will happen next
  all_results = matrix(0,nrow = n_sim * max_n_per_dec, ncol = 6);
  colnames(all_results) = c("sim_id","subj_id","dose_num","tox_prob","tox","stage");

  all_results[,"sim_id"] = rep(sim_id, each = max_n_per_dec);
  all_results[,"subj_id"] =
    unlist(sapply(sim_specific_start_id - 1,rep,each = max_n_per_dec, simplify = F)) +
    unlist(sapply(diff(c(0,which(!duplicated(all_results[,"sim_id"], fromLast=T)))), seq, from = 1,simplify = F))

  all_results[,"dose_num"] = rep(sim_specific_dose_start, each = max_n_per_dec);
  all_results[,"tox_prob"] = true_tox_curve[all_results[,"dose_num"]];
  all_results[,"tox"] = rbinom(nrow(all_results),1,all_results[,"tox_prob"]);
  all_results[,"stage"] = stage_label;


  #Now implement de-escalation rules, which
  #are assumed to act locally on each DEC or globally over
  #all DECs.
  remove_rows = NULL;
  if(module_rule == "local") {

    #SimID-specific list of numbers of subjects to look over
    #Initially equal to n
    num_at_final_level = rep(max_n_per_dec, n_sim);
    names(num_at_final_level) = sim_id;
    #The code uses `index' to look down each DEC
    #and determine when dose would have been de-escalated.
    index = seq_len(nrow(all_results));
    #Calculate cumulative toxicity rates after each patient
    tox_sums = tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),cumsum);
    tox_rates = lapply(tox_sums,"/",1:max_n_per_dec);
    #Before patient 'first_patient_look', we will only if stop if the absolute number of toxicities
    #exceeds the threshold rate *at* 'first_patient_look', i.e. if > floor(local_thresh * first_patient_look)
    if(first_patient_look > 1) {
      tox_rates = lapply(tox_rates,"[",-(1:(first_patient_look)));
      tox_sums = lapply(tox_sums,"[",(1:first_patient_look));
    }
    #If rate exceeds 'local_thresh', will de-escalate at the next patient.
    #This equals Inf if there is never a need to de-escalate
    decrease_after = pmin(first_patient_look + unlist(lapply(lapply(lapply(tox_rates,">",thresh_decrease),which),min,Inf)),
                          unlist(lapply(lapply(lapply(tox_sums,">",floor(thresh_decrease * first_patient_look)),which),min,Inf)));
    #Keep looping through until no de-escalation occurs
    while(any(decrease_after < num_at_final_level)) {
      #Sims to look at
      sim_set_char = names(decrease_after)[which(decrease_after < num_at_final_level)];
      sim_set = as.numeric(sim_set_char);
      #Determine whether this is the first run through the loop
      if(all(num_at_final_level == max_n_per_dec)) {
        #Start will all sim ids that need de-escalating
        index = which(all_results[,"sim_id"] %in% sim_set);
      } else {
        #Otherwise prune away the sim ids that have been handled by an earlier run
        index = setdiff(index,which(!all_results[,"sim_id"] %in% sim_set));
      }
      for(s in sim_set) {
        #For a given sim id, decrease index gives the set of rows,
        #the dose assignments of which would have been de-escalated
        #before the subjects are enrolled. It will be as if their
        #toxicity outcomes at the un-de-escalated dose were never observed
        #keep index gives the complement, the rows that are left alone
        decrease_index = intersect(index,
                                   which(all_results[,"sim_id"] == s));
        keep_index = decrease_index[which(all_results[decrease_index,"subj_id"] <= all_results[decrease_index,"subj_id"][decrease_after[as.character(s)]])];
        decrease_index = decrease_index[which(all_results[decrease_index,"subj_id"] > all_results[decrease_index,"subj_id"][decrease_after[as.character(s)]])];
        #if(sum(duplicated(all_results[decrease_index,"sim_specific_dose_start"]))!=length(decrease_index)-1) {stop();}
        #De-escalate
        stopifnot(length(decrease_index) || length(keep_index))
        all_results[decrease_index,"dose_num"] = all_results[decrease_index,"dose_num"] - 1;
        if(all_results[decrease_index,"dose_num"][1] > 0) {
          #assign new toxicity prob
          all_results[decrease_index,"tox_prob"] = true_tox_curve[all_results[decrease_index,"dose_num"]];
          #Determine toxicities under the de-escalated dose.
          all_results[decrease_index,"tox"] = rbinom(length(decrease_index),1,all_results[decrease_index,"tox_prob"]);
          #Prune away the keep index rows from the global index.
          index = setdiff(index,keep_index);
        } else {
          #Keep track of rows that would never have occurred because
          #de-escalation went down to zero
          remove_rows = c(remove_rows,decrease_index);
          #Remove this entire sim id (for this given dec id) from
          #consideration, as it has been addressed
          index = setdiff(index,which(all_results[,"sim_id"] == s));
          decrease_after[as.character(s)] = Inf;
        }
      }
      sim_set_char = setdiff(sim_set_char, unique(all_results[remove_rows,"sim_id"]));
      #Recalculate the list length to reflect the patients/sims that have been addresed
      #This should be monotonically decreasing for each run of the while loop
      num_at_final_level[sim_set_char] = tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),length)
      #Recalculate the toxicity rates
      tox_sums = tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),cumsum);
      tox_rates = mapply("/",tox_sums,sapply(num_at_final_level[sim_set_char],seq,from=1,simplify=F),SIMPLIFY=F);
      #Before patient 'first_patient_look', we will only if stop if the absolute number of toxicities
      #exceeds the threshold rate *at* 'first_patient_look', i.e. if > floor(local_thresh * first_patient_look)
      if(first_patient_look > 1) {
        tox_rates = lapply(tox_rates,"[",-(1:(first_patient_look)));
        tox_sums = lapply(tox_sums,"[",(1:first_patient_look));
      }
      #recalculate when to de-escalate the next time, if applicable
      decrease_after[sim_set_char] = pmin(first_patient_look + unlist(lapply(lapply(lapply(tox_rates,">",thresh_decrease),which),min,Inf)),
                                          unlist(lapply(lapply(lapply(tox_sums,">",floor(thresh_decrease * first_patient_look)),which),min,Inf)));
    }
  }
  #Remove those rows that would never have been enrolled because the
  #DEC(s) stopped for toxicity
  if(length(remove_rows)) {
    all_results = all_results[-remove_rows,,drop=F];
  }

  #Renumber the subjects in 'all_results' to reflect that the patients
  #in the `remove rows' set would never have been enrolled in the first place
  #Subtle note: it is crucial that the 'stage' indicator
  #be set before the subject ids are relabeled:
  skip_id = 1 + which(diff(all_results[,"subj_id"])>1);
  while(length(skip_id>0)) {
    all_results[skip_id,"subj_id"] = all_results[skip_id,"subj_id"] - 1;
    skip_id = 1 + which(diff(all_results[,"subj_id"])>1);
  }


  last_subj_data = all_results[!duplicated(all_results[,"sim_id"],fromLast = T),,drop = F];
  estMTD[sim_id] = ifelse(last_subj_data[,"subj_id"] < sim_specific_start_id + max_n_per_dec - 1, 0,
                          ifelse(decrease_after < Inf, last_subj_data[,"dose_num"] - 1, last_subj_data[,"dose_num"]));
  enrollment[sim_id] = tapply(all_results[,"subj_id"],list(all_results[,"sim_id"]),length)

  list(all_results = all_results[order(all_results[,"sim_id"],all_results[,"subj_id"]),],
       enrollment = enrollment,
       estMTD = estMTD,
       seed = seed);
}
