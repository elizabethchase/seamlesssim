#' This function simulates a specified number of seamless trials for each design configuration
#' provided.
#'
#' twostage_simulator is the primary workhorse of seamlesssim. It simulates complex seamless Phase I/II
#' oncology trials as discussed in the article by Boonstra et al. (Arxiv, 2020). It allows clinical
#' trialists to determine operating characteristics of trials that assess both
#' toxicity and efficacy with a range of different design and analytic approaches. For more detailed
#' information, see Boonstra et al. (Arxiv, 2020) and the vignette.
#'
#' @param array_id A positive integer identifier that will be appended as a column, without
#' modification, to all results. This is meant to be helpful to the user when calling this
#' function multiple times, e.g. in parallel.
#' @param n_sim A positive integer indicating how many simulated trials to conduct for each
#' design configuration.
#' @param primary_objectives A list containing three named elements: tox_target, tox_delta_no_exceed,
#' and eff_target, such that tox_target is between 0 and 1, tox_delta_no_exceed is between 0 and
#' (1 - tox_target), and eff_target is between 0 and 1. These choices delineate the primary
#' objectives of all designs to be simulated.  The true MTD is defined as the dose level with
#' true probability of DLT closest to tox_target but not exceeding tox_target + tox_delta_no_exceed,
#' and true acceptable dose level(s) are defined as any dose level that is less than or equal to
#' the MTD with a true probability of response at least as large as eff_target. Each design will
#' recommend the estimated MTD if its estimated efficacy probability is at least eff_target,
#' and otherwise recommend no dose level.
#' @param dose_outcome_curves  A list containing three named elements and an optional fourth
#' element: tox_curve, eff_curve, scenario, and, optionally, eff_curve_stage2. tox_curve is
#' the true toxicity curve of the doses; eff_curve is the true efficacy curve of the doses;
#' scenario is an identifier of which true data-generating scenario is being run (meant to be
#' helpful to the user when calling this function multiple times for different data-generating
#' scenarios).
#' @param design_list A list specifying all specific module choices. It will be a list of lists
#' of lists. The highest level of the list corresponds to each overall design to to be
#' evaluated; this should be as long as the number of designs that the user wants to compare.
#' The next level of the list gives the list of module choices for each design. It must have a named
#' component module1 and will optionally have named components module2...module4, taken from
#' the bolded values of Figure 1 in the manuscript referenced above. If any of module2 to
#' module4 are not provided, they are assumed to correspond to a choice of module2 =
#' list(name = "none"). Finally, the lowest level of the list gives the list of choices for each particular
#' module. Each list must have one entry named "name" to indicate the choice of module, and
#' also a value for every argument that is specific to that module. See the vignette for
#' examples.
#' @param stan_args A list containing eight named elements for the Bayesian isotonic regression.
#' For users without familiarity with STAN, stan_args can be left as NA (the default), and the defaults will all
#' be used. Alternatively, users can modify any/all of these arguments, leaving the others as defaults
#' or NA:
#' \describe{
#'    \item{n_mc_warmup}{A positive integer giving the number of desired warmup runs; the default is 1000}
#'    \item{n_mc_samps}{A positive integer giving the number of additional samples to run
#'    after warmup is completed; the default is 2000}
#'    \item{mc_chains}{A positive integer indicating the number of chains to run in parallel,
#'    which will multiply the final number of samples; the default is 4}
#'    \item{mc_thins}{A positive integer indicating the number of iterations to thin by
#'    (increasing thinning will decrease the final number of samples); the default is 1}
#'    \item{mc_stepsize}{A numeric value between 0 and 1 that is passed to control
#'    in the call to stan() as the stepsize argument; the default is 0.1}
#'    \item{mc_adapt_delta}{A numeric value between 0 and 1 that is passed to control
#'    in the call to stan() as the adapt_delta argument; the default is 0.8}
#'    \item{mc_max_treedepth}{A positive integer passed to control in the call to stan()
#'    as the max_treedepth argument; the default is 15}
#'    \item{ntries}{A positive integer. The stan algorithm throws warnings about divergent
#'    transitions, which are indicative of an inability to fully explore the posterior space.
#'    Sometimes this number can be extremely large, which suggests that the fitted model needs
#'    to be reparametrized. However, in this case, divergent transitions seem to be sporadic.
#'    ntries indicates how many reruns of the algorithm should be tried when > 0 divergent
#'    transitions are encountered. The run with the fewest such transitions is kept. The default is 2.}
#' }
#' For users without familiarity with STAN who still wish to use Bayesian isotonic regression,
#' some or all of these arguments may be left as NA, and default specifications will be used.
#' @param sim_labels A vector of anything but must be as long as n_sim. It will be included
#' in the final data.frame of results under a column name of sim_id. It is provided to allow
#' the user to uniquely identify simulations and is useful when this function is used in parallel.
#' @param design_labels A vector of anything but must be as long as length(design_list). It
#' will be included in the final data.frame of results under a column name of design. It is
#' provided to allow the user to uniquely identify designs and is useful when this function
#' is used in parallel.
#' @param do_efficient_simulation If TRUE, the simulator will run in such a way that, to the
#' maximum possible extent, simulated data will be reused between consecutive designs. So, for
#' example, design 1 may be identical to design 2 up to module 3, in which case the data
#' from modules 1 and 2 can be reused from design 1 to design 2. If FALSE, each design will be
#' simulated independently of each other design, but the whole simulator will take longer to run.
#' @param random_seed A positive integer seed set prior to starting the simulations.
#' @return The function returns a named list with entries:
#' \describe{
#'    \item{patient_data}{A data.frame with number of rows equal to number of individual patients
#'    simulated across all simulations of all designs, i.e. if every single design were to enroll
#'    the maximum possible number of patients, say, n, the number of rows would be n *
#'    length(design_list) * n_sim. }
#'    \item{sim_data_stage1}{A data.frame with number of rows equal to length(design_list) * n_sim,
#'     i.e. one per design per simulation. It gives trial-level summary information about the
#'     status of the trial at the end of module 2 of each design.}
#'    \item{sim_data_stage2}{A data.frame with number of rows equal to length(design_list) * n_sim,
#'     i.e. one per design per simulation. It gives trial-level summary information about the
#'     status of the trial at the end of module 4 of each design.}
#'     \item{dose_outcome_curves}{The user-inputted argument to this function having the same name.}
#'     \item{titecrm_args}{The list of common arguments that were used for the crm simulator.}
#'     \item{design_list}{The user-inputted argument to this function having the same name.}
#'     \item{design_description}{A character matrix with number of rows equal to length(design_list)
#'      and number of columns equal to the total number of modules used in the trial, presumably
#'      4. It is meant to give a concise, simple summary and comparison of each design, without
#'      going into the details of each design.}
#'      \item{shared_design_elements}{An integer matrix with number of rows equal to length(design_list)
#'      and number of columns equal to the total number of modules used in the trial, presumably
#'      4. It gives the simulators assessment of which design elements could be recycled (therefore
#'      saving time if do_efficient_simulation==TRUE).}
#'      \item{random_seed}{The user-inputted argument to this function having the same name.}
#' }
#' @importFrom dplyr summarise near %>% bind_cols arrange group_by summarize mutate select pull
#' left_join
#' @importFrom stats rbinom qbeta xtabs
#' @import binom
#' @import rstan
#' @export
twostage_simulator = function(array_id = 1,
                              n_sim,
                              primary_objectives,
                              dose_outcome_curves,
                              design_list,
                              stan_args = NA,
                              sim_labels = NULL,
                              design_labels = NULL,
                              do_efficient_simulation = T,
                              verbose = F,
                              random_seed = 1) {

  # Setup----
  y <- NULL
  sim_id <- NULL
  subj_id <- NULL
  stage <- NULL
  dose_num <- NULL
  estMTD <- NULL
  receivedEstMTD <- NULL
  eff <- NULL
  model_mean_prob <- NULL
  model_lower_ci_prob <- NULL
  scenario <- NULL
  design <- NULL

  eps = .Machine$double.eps^0.75;
  n_dose = length(dose_outcome_curves[["tox_curve"]]);
  if(is.null(sim_labels)) {sim_labels = 1:n_sim;}
  stopifnot(near(length(sim_labels), n_sim));
  stopifnot("list" %in% class(design_list));
  if(is.null(design_labels)) {design_labels = 1:length(design_list);}
  stopifnot(near(length(design_labels), length(design_list)));
  if(primary_objectives[["tox_target"]] > 1 | primary_objectives[["tox_target"]] < 0){
    warning("tox_target in primary_objectives is outside of [0,1]")
  }
  if(primary_objectives[["tox_delta_no_exceed"]] > (1-primary_objectives[["tox_target"]]) | primary_objectives[["tox_delta_no_exceed"]] < 0){
    warning("tox_delta_no_exceed in primary_objectives is outside of [0,1-tox_target]")
  }
  if(primary_objectives[["eff_target"]] > 1 | primary_objectives[["eff_target"]] < 0){
    warning("eff_target in primary_objectives is outside of [0,1]")
  }

  # + Dose-Efficacy curves----
  #Format true dose-efficacy curves, which are allowed to differ [in truth] between stage 1 ("escalation") and stage 2 ("expansion");
  #Also recorded here (i) are the doses that are 'acceptable' (exceeding eff_target in efficacy, and falling below tox_target + tox_delta_no_exceed) in toxicity
  #and (ii) the 'preferred' dose (the largest acceptable dose)
  store_eff_curves =
    matrix(0, 2, (2 * n_dose) + 1,
           dimnames = list(c("stage1","stage2"),
                           c(paste0("dose",1:n_dose),paste0("dose",1:n_dose,"_accept"),"pref_dose")));
  if(!all(c("tox_curve","eff_curve","scenario") %in% names(dose_outcome_curves))) {
    stop("'dose_outcome_curves' must be a named list containing named entries 'tox_curve', 'eff_curve', 'scenario'");
  }
  if(!near(length(dose_outcome_curves[["tox_curve"]]),length(dose_outcome_curves[["eff_curve"]]))) {
    stop("The vectors 'tox_curve' and 'eff_curve', which are elements of 'dose_outcome_curves', must have the same length");
  }
  if(any(dose_outcome_curves[["tox_curve"]] > 1) | any(dose_outcome_curves[["tox_curve"]] < 0)){
    warning("tox_curve has probabilities outside of [0, 1]")
  }
  if(any(dose_outcome_curves[["eff_curve"]] > 1) | any(dose_outcome_curves[["eff_curve"]] < 0)){
    warning("eff_curve has probabilities outside of [0, 1]")
  }

  for(curr_stage in c("stage1","stage2")) {
    #The same true efficacy curve applies to each stage of the trial unless a named element 'eff_curve_stage2'
    #says otherwise
    if(curr_stage == "stage1") {
      curr_efficacy_curve = dose_outcome_curves[["eff_curve"]];
    } else if("eff_curve_stage2" %in% names(dose_outcome_curves)) {
      if(!near(length(dose_outcome_curves[["eff_curve"]]),length(dose_outcome_curves[["eff_curve_stage2"]]))) {
        stop("The vectors 'eff_curve' and 'eff_curve_stage2', which are elements of 'dose_outcome_curves', must have the same length");
      }
      if(any(dose_outcome_curves[["eff_curve_stage2"]] > 1) | any(dose_outcome_curves[["eff_curve_stage2"]] < 0)){
        warning("eff_curve_stage2 has probabilities outside of [0, 1]")
      }
      curr_efficacy_curve = dose_outcome_curves[["eff_curve_stage2"]];
    }

    #The first 'T' corresponds to dose zero, which is acceptable if andonly if all other dose levels are unacceptable
    curr_acceptability =
      c(T,(curr_efficacy_curve >= (primary_objectives["eff_target"])) &
          (dose_outcome_curves[["tox_curve"]] <= primary_objectives["tox_target"] + primary_objectives["tox_delta_no_exceed"]));
    curr_pref = max(which(curr_acceptability)) - 1;
    store_eff_curves[curr_stage,] = c(curr_efficacy_curve,curr_acceptability[-1],curr_pref);
    rm(curr_acceptability, curr_pref);
  }
  rm(curr_efficacy_curve, curr_stage);

  #Calculate the true MTD, which is constant for the entire simulation
  trueMTD = max(c(which.min(abs(dose_outcome_curves[["tox_curve"]] - primary_objectives["tox_target"]) /
                              (dose_outcome_curves[["tox_curve"]] <= (primary_objectives["tox_target"] + primary_objectives["tox_delta_no_exceed"]))),0));

  # + Error Checking----
  #Concise description matrix of each design using the 'name' component of each module
  ncol_needed =  4
  design_description = matrix(NA, nrow = length(design_list),
                              ncol = ncol_needed,
                              dimnames = list(NULL, c("module1","module2","module3","module4")));

  # Check for coherent combinations of module choices.
  for(i in seq_along(design_list)) {
    # ++ Module names----
    #Module1 must be specified by user; other modules may be skipped. Modules that are not provided
    #are assumed to be not wanted and set to 'none'.
    foo = which(!c("module1","module2","module3","module4") %in% names(design_list[[i]]));
    for(j in 1:4) {
      #This takes care of the modules that were not provided at all
      if(j %in% foo) {
        if(near(j, 1)) {
          stop(paste0("Error in entry ",i," of 'design_list': 'module1' was not found but must be a named list"));
        } else {
          design_list[[i]] = c(design_list[[i]],list(list("name" = "none")));
          names(design_list[[i]])[length(design_list[[i]])] = paste0("module",j);
        }
      }
      #This takes care of the modules that were provided as empty lists
      if(isTRUE(all.equal(list(),design_list[[i]][[paste0("module",j)]]))) {
        if(near(j, 1)) {
          stop(paste0("Error in entry ",i," of 'design_list': 'module1' was not found but must be a named list"));
        } else {
          design_list[[i]][[paste0("module",j)]] = list("name" = "none")
        }
      }
      if(!"name" %in% names(design_list[[i]][[paste0("module",j)]])) {
        stop(paste0("Error in entry ",i," of 'design_list': no name was provided for 'module",j,"' method"));
      }
      design_list[[i]] = design_list[[i]][order(names(design_list[[i]]))];
    }
    #
    design_description[i,c("module1","module2","module3","module4")] =
      c(design_list[[i]][["module1"]]$name,
        design_list[[i]][["module2"]]$name,
        design_list[[i]][["module3"]]$name,
        design_list[[i]][["module4"]]$name);
    #Check for valid module1 name
    if(!design_description[i,"module1"] %in% c("crm", "3pl3", "empiric","fixed")) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module1' is '", design_description[i,"module1"] ,"' but must be 'crm', '3pl3', 'empiric', or 'fixed'"));
    }
    #Check for valid module2 name
    if(!design_description[i,"module2"] %in% c("none", "bayes", "bayes_isoreg", "inverted_score", "min_num_resp", "min_pct_resp")) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' is '", design_description[i,"module2"], "' but  must be 'none', 'bayes', 'bayes_isoreg', 'inverted_score', 'min_num_resp', or 'min_pct_resp'"));
    }
    #Check for valid module3 name
    if(!design_description[i,"module3"] %in% c("crm", "continue_crm", "3pl3", "empiric", "fixed","none")) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module3' is '", design_description[i,"module3"] ,"' but must be 'crm', 'continue_crm', '3pl3', 'empiric', 'fixed', or 'none'"));
    }
    #Check for valid module4 name
    if(!design_description[i,"module4"] %in% c("none", "bayes", "bayes_isoreg", "inverted_score", "min_num_resp", "min_pct_resp")) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' is '", design_description[i,"module4"], "' but  must be 'none', 'bayes', 'bayes_isoreg', 'inverted_score', 'min_num_resp', or 'min_pct_resp'"));
    }

    #Check for valid module1+module3 combination
    if(design_description[i,"module1"] != "crm" &&
       design_description[i,"module3"] %in% "continue_crm") {
      stop(paste0("Error in entry ",i," of 'design_list': 'module1' is '", design_description[i,"module1"] ,"' but must be 'crm' when 'module3' is 'continue_crm'"));
    }
    #Check for valid module1+module2 combination
    if(design_description[i,"module1"] == "3pl3" &&
       design_description[i, "module2"] == "min_num_resp"){
      warning(paste0("In entry ",i," of 'design_list': using '3pl3' for module1 followed by 'min_num_resp' for module2 is not advised"))
    }
    #Check for valid module3+module4 combination
    if(design_description[i,"module3"] == "3pl3" &&
       design_description[i, "module4"] == "min_num_resp"){
      warning(paste0("In entry ",i," of 'design_list': using '3pl3' for module3 followed by 'min_num_resp' for module4 is not advised"))
    }
    # ++ Module 1----
    if(design_description[i,"module1"] == "crm" &&
       length(foo <- setdiff(c("n","skeleton","starting_dose","beta_scale","dose_cohort_size","dose_cohort_size_first_only","earliest_stop"),
                             names(design_list[[i]][["module1"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module1' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module1"] == "3pl3" &&
       length(foo <- setdiff(c("starting_dose"),
                             names(design_list[[i]][["module1"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module1' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module1"] == "3pl3" &&
       length(foo <- setdiff(names(design_list[[i]][["module1"]]),
                             c("starting_dose", "name")))) {
      warning(paste0("Entry ",i," of 'design_list': 'module1' has unnecessary components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module1"] == "empiric" &&
       length(foo <- setdiff(c("n","starting_dose","rule","first_patient_look","thresh_decrease"),
                             names(design_list[[i]][["module1"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module1' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module1"] == "fixed" &&
       length(foo <- setdiff(c("n","starting_dose"),
                             names(design_list[[i]][["module1"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module1' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module1"] == "fixed" &&
       length(foo <- setdiff(names(design_list[[i]][["module1"]]),
                             c("n", "starting_dose", "name")))) {
      warning(paste0("Entry ",i," of 'design_list': 'module1' has unnecessary components: ", paste0(foo,collapse=", ")));
    }

    # ++ Module 2----
    if(design_description[i,"module2"] == "bayes" &&
       length(foo <- setdiff(c("prob_threshold","prior_mean","prior_n_per"),
                             names(design_list[[i]][["module2"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module2"] == "bayes" &&
       !near(length(design_list[[i]][["module2"]]$prior_mean),
             length(dose_outcome_curves[["tox_curve"]]))) {
      stop(paste0("Error in entry ",i," of 'design_list': in 'module2', the length of 'prior_mean' is not equal to the implied number of dose levels, which is ",length(dose_outcome_curves[["tox_curve"]])));
    }
    if(design_description[i,"module2"] == "bayes" &&
       (design_list[[i]][["module2"]]$prior_n_per <= 0)) {
      stop(paste0("Error in entry ",i," of 'design_list': in 'module2', 'prior_n_per' is non-positive but must be positive"));
    }
    if(design_description[i,"module2"] == "bayes" &&
       (any(design_list[[i]][["module2"]]$prior_mean < 0) ||
        any(design_list[[i]][["module2"]]$prior_mean > 1))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' has an element of 'prior_mean' outside of the interval (0,1), which is not allowed"));
    }
    if(design_description[i,"module2"] == "bayes_isoreg" &&
       length(foo <- setdiff(c("prob_threshold","alpha_scale"),
                             names(design_list[[i]][["module2"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module4"] == "bayes_isoreg" &&
       design_list[[i]][["module4"]]$alpha_scale <= 0) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' must have 'alpha_scale' > 0"));
    }
    if(design_description[i,"module2"] %in% "inverted_score" &&
       length(foo <- setdiff(c("ci_level_onesided"),
                             names(design_list[[i]][["module2"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module2"] %in% "inverted_score" &&
       (design_list[[i]][["module2"]]$ci_level_onesided < 0.5 ||
        design_list[[i]][["module2"]]$ci_level_onesided > 1.0)) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' has 'ci_level_onesided' outside of the interval [0.5, 1], which is not allowed"));
    }
    if(design_description[i,"module2"] == "min_num_resp" &&
       !"number"%in%names(design_list[[i]][["module2"]])) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' is missing the named component 'number'"));
    }
    if(design_description[i,"module2"] == "min_num_resp" &&
       design_list[[i]][["module2"]]$number > design_list[[i]][["module1"]]$n) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' would require more responders than maximum planned stage 1 enrollment"));
    }
    if(design_description[i,"module2"] == "min_pct_resp" &&
       !"percent"%in%names(design_list[[i]][["module2"]])) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module2' is missing the named component 'percent'"));
    }
    # ++ Module 3----
    if(design_description[i,"module3"] == "crm") {
    #all crm ingredients must be provided otherwise
      foo <- setdiff(c("n","skeleton","beta_scale","dose_cohort_size","dose_cohort_size_first_only","earliest_stop"),
                       names(design_list[[i]][["module3"]]))

      if(length(foo)) {
        stop(paste0("Error in entry ",i," of 'design_list': 'module3' is missing the following required components: ", paste0(foo,collapse=", ")));
      }
    }
    if(design_description[i,"module3"] == "continue_crm" &&
       length(foo <- setdiff("n",
                             names(design_list[[i]][["module3"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module3' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module3"] == "3pl3" &&
       length(setdiff(foo <- names(design_list[[i]][["module3"]]),
                      "name"))) {
      warning(paste0("Entry ",i," of 'design_list': 'module3' has unnecessary components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module3"] == "empiric" &&
       length(foo <- setdiff(c("n","rule","first_patient_look","thresh_decrease"),
                             names(design_list[[i]][["module3"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module3' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module3"] == "fixed" &&
       length(foo <- setdiff("n",
                             names(design_list[[i]][["module3"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module3' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module3"] == "fixed" &&
       length(foo <- setdiff(names(design_list[[i]][["module3"]]),
                             c("n", "name")))) {
      warning(paste0("Note that in entry ",i," of 'design_list': 'module3' has unnecessary components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module3"] == "none") {
      design_list[[i]][["module3"]]$n = 0;
    }
    # ++ Module 4----
    if(design_description[i,"module4"] == "bayes" &&
       length(foo <- setdiff(c("prob_threshold","prior_mean","prior_n_per","include_stage1_data"),
                             names(design_list[[i]][["module4"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module4"] == "bayes" &&
       !near(length(design_list[[i]][["module4"]]$prior_mean),
             length(dose_outcome_curves[["tox_curve"]]))) {
      stop(paste0("Error in entry ",i," of 'design_list': in 'module4', the length of 'prior_mean' is not equal to the implied number of dose levels, which is ",length(dose_outcome_curves[["tox_curve"]])));
    }
    if(design_description[i,"module4"] == "bayes" &&
       (design_list[[i]][["module4"]]$prior_n_per <= 0)) {
      stop(paste0("Error in entry ",i," of 'design_list': in 'module4', 'prior_n_per' is non-positive but must be positive"));
    }
    if(design_description[i,"module4"] == "bayes" &&
       (any(design_list[[i]][["module4"]]$prior_mean < 0) ||
        any(design_list[[i]][["module4"]]$prior_mean > 1))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' has an element of 'prior_mean' outside of the interval (0,1), which is not allowed"));
    }
    if(design_description[i,"module4"] == "bayes_isoreg" &&
       length(foo <- setdiff(c("prob_threshold","alpha_scale","include_stage1_data"),
                             names(design_list[[i]][["module4"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module4"] == "bayes_isoreg" &&
       design_list[[i]][["module4"]]$alpha_scale <= 0) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' must have 'alpha_scale' > 0"));
    }
    if(design_description[i,"module4"] == "inverted_score" &&
       length(foo <- setdiff(c("ci_level_onesided","include_stage1_data"),
                             names(design_list[[i]][["module4"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module4"] %in% "inverted_score" &&
       (design_list[[i]][["module4"]]$ci_level_onesided < 0.5 ||
        design_list[[i]][["module4"]]$ci_level_onesided > 1.0)) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' has 'ci_level_onesided' outside of the interval [0.5, 1], which is not allowed"));
    }
    if(design_description[i,"module4"] == "min_num_resp" &&
       length(foo <- setdiff(c("number","include_stage1_data"),
                             names(design_list[[i]][["module4"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module4"] == "min_num_resp" &&
       design_list[[i]][["module4"]]$number > (design_list[[i]][["module4"]]$include_stage1_data * design_list[[i]][["module1"]]$n) + design_list[[i]][["module4"]]$n) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' would require more responders than maximum planned total enrollment"));
    }
    if(design_description[i,"module4"] == "min_pct_resp" &&
       length(foo <- setdiff(c("percent","include_stage1_data"),
                             names(design_list[[i]][["module4"]])))) {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' is missing the following required components: ", paste0(foo,collapse=", ")));
    }
    if(design_description[i,"module4"] != "none" &&
       !design_list[[i]][["module4"]]$include_stage1_data &&
       design_description[i,"module3"] == "none") {
      stop(paste0("Error in entry ",i," of 'design_list': 'module4' must have 'include_stage1_data = T' since no second stage data is generated, i.e. module3 is not run"));
    }

  }
  rm(i);

  # + Reordering of designs----
  design_description = data.frame(design_description);
  design_description[,"module1"] = factor(design_description[,"module1"], levels = c("crm","3pl3","empiric","fixed"), ordered = T);
  design_description[,"module2"] = factor(design_description[,"module2"], levels = c( "bayes", "bayes_isoreg", "inverted_score", "min_num_resp","min_pct_resp","none"), ordered = T)
  design_description[,"module3"] = factor(design_description[,"module3"], levels = c("continue_crm", "crm", "3pl3", "empiric", "fixed", "none"), ordered = T)
  design_description[,"module4"] = factor(design_description[,"module4"], levels = c("bayes", "bayes_isoreg", "inverted_score", "min_num_resp","min_pct_resp", "none"), ordered = T)
  if(any(is.na(design_description))) {
    which_designs <- which(rowSums(is.na(design_description)) > 0);
    stop(paste0("There was a problem: one of the modules in design(s), i.e. ",paste0(which_designs,collapse = " and "), ", has an unrecognized name"));
    rm(which_designs);
  }

  # If efficient simulations are requested, internally reorder the provided design lists so as to most efficiently re-use the
  # simulated data between designs.
  if(do_efficient_simulation) {
    design_list_reorder =
      order(design_description[,"module1"],
            design_description[,"module3"],
            design_description[,"module2"],
            design_description[,"module4"]);
  } else {
    design_list_reorder = seq_len(nrow(design_description));
  }

  design_list = design_list[design_list_reorder];
  design_description = design_description[design_list_reorder,];

  shared_design_elements = matrix(seq_along(design_list),
                                  nrow = length(design_list),
                                  ncol = 4,
                                  dimnames = list(NULL,paste0("module",c(1:4))));
  # This identifies common elements between designs, allowing the simulator
  # to reuse simulated data between designs and thereby be more efficient.
  if(do_efficient_simulation) {
    for(i in seq_along(design_list)[-1]) {
      for(j in seq_len(i-1)) {
        curr_match = shared_design_elements[i,];
        for(k in colnames(shared_design_elements)) {
          if(isTRUE(all.equal(design_list[[i]][[k]],
                              design_list[[j]][[k]]))) {
            curr_match[k] = shared_design_elements[j,k];
          }
          # crm -> none is equivalent to crm -> continue_crm
          # (with regard to toxicity outcomes) when total sample size is
          # identical.
          if(k == "module3" &&
             design_list[[i]][["module1"]]$name == "crm" &&
             design_list[[j]][["module1"]]$name == "crm" &&
             isTRUE(all.equal(design_list[[i]][["module1"]][c("skeleton","starting_dose","beta_scale","dose_cohort_size","dose_cohort_size_first_only","earliest_stop")],
                              design_list[[j]][["module1"]][c("skeleton","starting_dose","beta_scale","dose_cohort_size","dose_cohort_size_first_only","earliest_stop")])) &&
             ((design_list[[i]][["module3"]]$name == "none" &&
               design_list[[j]][["module3"]]$name == "none" &&
               near(design_list[[i]][["module1"]]$n,
                    design_list[[j]][["module1"]]$n)) ||
              (design_list[[i]][["module3"]]$name == "none" &&
               design_list[[j]][["module3"]]$name == "continue_crm" &&
               near(design_list[[i]][["module1"]]$n,
                    design_list[[j]][["module1"]]$n + design_list[[j]][["module3"]]$n)) ||
              (design_list[[i]][["module3"]]$name == "continue_crm" &&
               design_list[[j]][["module3"]]$name == "continue_crm" &&
               near(design_list[[i]][["module1"]]$n + design_list[[i]][["module3"]]$n,
                    design_list[[j]][["module1"]]$n + design_list[[j]][["module3"]]$n)))) {
            curr_match["module1"] = shared_design_elements[j,"module1"];
            curr_match["module3"] = shared_design_elements[j,"module3"];
          }
        }
        if(sum(curr_match < i) > sum(shared_design_elements[i,] < i)) {
          shared_design_elements[i,] = curr_match;
        }
      }
      rm(curr_match)
    }
    rm(i,j, k);
  }

  if(any(design_description[,"module2"] == "bayes_isoreg") ||
     any(design_description[,"module4"] == "bayes_isoreg")) {
    #Fill in STAN default parameters

    if (all(is.na(stan_args))) {
      stan_args <- list(
        n_mc_warmup = 1e3,
        n_mc_samps = 2e3,
        mc_chains = 4,
        mc_thin = 1,
        mc_stepsize = 0.1,
        mc_adapt_delta = 0.8,
        mc_max_treedepth = 15,
        ntries = 2
      )
    } else if (is.na(stan_args$n_mc_warmup)){
      stan_args$n_mc_warmup <- 1e3
    } else if (is.na(stan_args$n_mc_samps)){
      stan_args$n_mc_samps <- 2e3
    } else if (is.na(stan_args$mc_chains)){
      stan_args$mc_chains <- 4
    } else if (is.na(stan_args$mc_thin)){
      stan_args$mc_thin <- 1
    } else if (is.na(stan_args$mc_stepsize)){
      stan_args$mc_stepsize <- 0.1
    } else if (is.na(stan_args$mc_adapt_delta)){
      stan_args$mc_adapt_delta <- 0.8
    } else if (is.na(stan_args$mc_max_treedepth)){
      stan_args$mc_max_treedepth <- 15
    } else if (is.na(stan_args$ntries)){
      stan_args$ntries <- 2
    }
  }

  store_module1_dat =
    store_stage1_estMTD =
    store_stage1_enrollment =
    store_stage1_summary =
    store_stage1_RP2D =
    # module3 data will only be observed if module2 is successfully passed
    store_module3_dat_premodule2 =
    # Thus we need two versions of the module3 data: pre and post module2. The
    # latter will be empty for trials that stop at module2.
    store_module3_dat_postmodule2 =
    store_stage2_estMTD_premodule2 =
    store_stage2_estMTD_postmodule2 =
    store_stage2_enrollment_premodule2 =
    store_stage2_enrollment_postmodule2 =
    vector("list", length(design_list));
  patient_data_to_return =
    sim_data_stage1 =
    sim_data_stage2 = NULL;

  # Write out the initial values of the 'titecrm_args', which will be passed to the 'titesim_phil' function.
  # Anything that is NA is as such because that is a design-specific choice that is filled in down below
  starting_titecrm_args = list(
    PI = dose_outcome_curves[["tox_curve"]],
    prior = NA,#this is 'skeleton' from each module
    target = primary_objectives[["tox_target"]],
    n = NA,#this is 'n' from each module
    x0 = NA,#this is the starting dose level and is either 'starting_dose' for module 1 or trial-dependent for module 3
    nsim = 1,#this is the number of simulations to do *within* the 'titesim_phil' function and so is always 1
    restrict = TRUE,#always enact dose-escalation constraints
    obswin = 1,#we're not simulating the tite-part of the trial, so this should stay fixed
    rate = 1/2,#we're not simulating the tite-part of the trial, so this should stay fixed
    accrual = "fixed",#we're not simulating the tite-part of the trial, so this should stay fixed
    count = FALSE,#don't gi
    method = "bayes",#only 'bayes' is implemented
    model = "empiric",#only 'empiric', i.e. the power model, is implemented
    scale = NA,#this is 'beta_scale' from reach module
    no.exceed = primary_objectives[["tox_delta_no_exceed"]],
    cohort.size = NA,#this is 'dose_cohort_size' from each module
    first.cohort.only = NA,#this is 'dose_cohort_size_first_only' from each module
    n.at.MTD = Inf#this is not used here
  );
  set.seed(random_seed);
  seeds_by_module = sample(.Machine$integer.max, 4);

  for(k in seq_along(design_list)) {
    #Start with a clean slate
    if("titecrm_args"%in%ls()) {rm(titecrm_args);}
    curr_design = design_list[[k]];

    #Module 1: Dose Assignments----
    set.seed(seeds_by_module[1]);
    #The following is only ever true if 'do_efficient_simulation == TRUE'. It
    # checks for shared elements across designs and resuses the already
    # simulated data as much as possible. Otherwise, each iteration completely
    # regenerates the data
    if(shared_design_elements[k,"module1"] != k) {
      if(verbose) {cat("reusing Module 1 conclusions...\n\n");}
      if(curr_design[["module1"]]$name == "3pl3") {
        module1_dat = store_module1_dat[[shared_design_elements[k,"module1"]]];
        module1_dat[,"design"] = design_list_reorder[k];
        stage1_estMTD = store_stage1_estMTD[[shared_design_elements[k,"module1"]]];
        stage1_enrollment = store_stage1_enrollment[[shared_design_elements[k,"module1"]]];
        # The below condition checks for the scenario that the identical dose assignment schemes have already been run
        # The first check is potentially subtly and warrants additional
        # explanation: the code has reordered the designs so that
        # crm (module1) -> continue_crm (module3) always comes before
        # crm (module1) -> none (module3). When the crm design settings are
        # equivalent, including total sample size, we can therefore reuse
        # the data from crm -> continue_crm for crm -> none. However, depending
        # on differences in module2, one design may actually stop at module2, and
        # so it would never see the module3 toxicity data.
      } else if(all(shared_design_elements[k,c("module1","module3")] != k) &&
                duplicated(shared_design_elements[1:k,paste0("module",c(1,3))])[k] &&
                design_list[[k]][["module1"]]$name == "crm" &&
                design_list[[k]][["module3"]]$name == "none" &&
                design_list[[shared_design_elements[k,"module1"]]][["module1"]]$name == "crm" &&
                design_list[[shared_design_elements[k,"module1"]]][["module3"]]$name == "continue_crm") {
        which_duplicate_of =
          which(rowSums(abs(shared_design_elements[1:(k-1),paste0("module",c(1,3)),drop = F] -
                              shared_design_elements[rep(k, k-1),paste0("module",c(1,3)),drop=F])) < eps)[1];
        module1_dat =
          rbind(store_module1_dat[[which_duplicate_of]],
                # "premodule2" means that these data may not end up being used
                # since the futility analysis may stop the trial at module2.
                store_module3_dat_premodule2[[which_duplicate_of]]) %>%
          arrange(sim_id,subj_id);
        module1_dat[,"stage"] = 1;
        module1_dat[,"design"] = design_list_reorder[k];
        stage1_estMTD = store_stage2_estMTD_premodule2[[which_duplicate_of]];
        stage1_enrollment = store_stage1_enrollment[[which_duplicate_of]] +
          store_stage2_enrollment_premodule2[[which_duplicate_of]];
        rm(which_duplicate_of);
      } else if(all(shared_design_elements[k,c("module1","module3")] != k) &&
                duplicated(shared_design_elements[1:k,paste0("module",c(1,3))])[k] &&
                design_list[[k]][["module1"]]$name == "crm" &&
                design_list[[k]][["module3"]]$name == "continue_crm" &&
                design_list[[shared_design_elements[k,"module1"]]][["module1"]]$name == "crm" &&
                design_list[[shared_design_elements[k,"module1"]]][["module3"]]$name == "continue_crm") {
        which_duplicate_of =
          which(rowSums(abs(shared_design_elements[1:(k-1),paste0("module",c(1,3)),drop = F] -
                              shared_design_elements[rep(k, k-1),paste0("module",c(1,3)),drop=F])) < eps)[1];
        module1_dat = store_module1_dat[[which_duplicate_of]];
        module1_dat[,"design"] = design_list_reorder[k];
        stage1_estMTD = store_stage1_estMTD[[which_duplicate_of]];
        stage1_enrollment = store_stage1_enrollment[[which_duplicate_of]];

        module3_dat = store_module3_dat_premodule2[[which_duplicate_of]];
        if(nrow(module3_dat) > 0) {
          module3_dat[,"design"] = design_list_reorder[k];
        }
        stage2_estMTD = store_stage2_estMTD_premodule2[[which_duplicate_of]];
        stage2_enrollment = store_stage2_enrollment_premodule2[[which_duplicate_of]];
        rm(which_duplicate_of);
      } else {
        module1_dat = store_module1_dat[[shared_design_elements[k,"module1"]]];
        module1_dat[,"design"] = design_list_reorder[k];
        stage1_estMTD = store_stage1_estMTD[[shared_design_elements[k,"module1"]]];
        stage1_enrollment = store_stage1_enrollment[[shared_design_elements[k,"module1"]]];
      }
    } else {
      if(curr_design[["module1"]]$name == "crm") {
        titecrm_args = starting_titecrm_args;
        titecrm_args$prior = curr_design[["module1"]]$skeleton;
        titecrm_args$scale = curr_design[["module1"]]$beta_scale;
        titecrm_args$x0 = curr_design[["module1"]]$starting_dose;
        titecrm_args$cohort.size = curr_design[["module1"]]$dose_cohort_size;
        titecrm_args$first.cohort.only = curr_design[["module1"]]$dose_cohort_size_first_only;
        titecrm_args$earliest_stop = curr_design[["module1"]]$earliest_stop;
        #Look ahead to see if entire CRM should be run all at once
        if(curr_design[["module3"]]$name == "continue_crm") {
          titecrm_args$n = curr_design[["module1"]]$n + curr_design[["module3"]]$n;
          second_stage_start_after = curr_design[["module1"]]$n;
        } else {
          titecrm_args$n = curr_design[["module1"]]$n;
          second_stage_start_after = Inf;
        }

        ##CRM toxicity data;
        foo = sim_twostage_crm(n_sim = n_sim,
                               titecrm_args = titecrm_args,
                               second_stage_start_after = second_stage_start_after);

        #The data at this point do not yet take into account stopping for futility, which is module-2-specific
        module1_dat = data.frame(foo[["all_results"]]);
        stage1_estMTD = foo[["first_stage_estMTD"]];
        stage1_enrollment = foo[["first_stage_enrollment"]];
        if(curr_design[["module3"]]$name == "continue_crm") {
          stage2_estMTD = foo[["second_stage_estMTD"]];
          stage2_enrollment = foo[["second_stage_enrollment"]];
          if(any(near(module1_dat[,"stage"], 2))) {
            module3_dat = cbind(array_id = array_id,
                                scenario = dose_outcome_curves[["scenario"]],
                                design = design_list_reorder[k],
                                filter(module1_dat, near(stage, 2)),
                                eff = NA,
                                eff_prob = NA);
            #Generate efficacy data
            module3_dat[,"eff_prob"] = store_eff_curves["stage2",module3_dat[,"dose_num"]];
            module3_dat[,"eff"] = rbinom(nrow(module3_dat),1,module3_dat[,"eff_prob"]);
          } else {
            module3_dat = data.frame(array_id = integer(0),
                                     scenario = integer(0),
                                     design = integer(0),
                                     filter(module1_dat, near(stage, 2)),
                                     eff = numeric(0),
                                     eff_prob = numeric(0));
          }
          store_module3_dat_premodule2[[k]] = module3_dat;
          store_stage2_estMTD_premodule2[[k]] = stage2_estMTD;
          store_stage2_enrollment_premodule2[[k]] = stage2_enrollment;
          module1_dat = filter(module1_dat, near(stage, 1));
        }
        rm(foo);
      } else if(curr_design[["module1"]]$name == "3pl3") {
        foo = sim_3pl3(n_sim = n_sim,
                       true_tox_curve = dose_outcome_curves[["tox_curve"]],
                       stage_label = 1,
                       sim_specific_start_id = NULL,
                       sim_specific_dose_start = rep(curr_design[["module1"]]$starting_dose, n_sim));

        module1_dat = data.frame(foo[["all_results"]]);
        stage1_estMTD = foo[["estMTD"]];
        stage1_enrollment = foo[["enrollment"]];
        rm(foo);
      } else if(curr_design[["module1"]]$name == "empiric") {

        foo = sim_empiric_dec(n_sim = n_sim,
                              true_tox_curve = dose_outcome_curves[["tox_curve"]],
                              stage_label = 1,
                              sim_specific_start_id = NULL,
                              sim_specific_dose_start = rep(curr_design[["module1"]]$starting_dose, n_sim),
                              max_n_per_dec = curr_design[["module1"]]$n,
                              module_rule = curr_design[["module1"]]$rule,
                              thresh_decrease = curr_design[["module1"]]$thresh_decrease,
                              first_patient_look = curr_design[["module1"]]$first_patient_look);
        module1_dat = data.frame(foo[["all_results"]]);
        stage1_estMTD = foo[["estMTD"]];
        stage1_enrollment = foo[["enrollment"]];
        rm(foo);
      } else if(curr_design[["module1"]]$name == "fixed") {
        foo = sim_empiric_dec(n_sim = n_sim,
                              true_tox_curve = dose_outcome_curves[["tox_curve"]],
                              stage_label = 1,
                              sim_specific_start_id = rep(1, n_sim),
                              sim_specific_dose_start = rep(curr_design[["module1"]]$starting_dose, n_sim),
                              max_n_per_dec = curr_design[["module1"]]$n,
                              module_rule = "local",
                              thresh_decrease = Inf)#Anything greater than or equal to 1 means never de-escalate;
        module1_dat = data.frame(foo[["all_results"]]);
        module1_dat[,"stage"] = 1;
        stage1_estMTD = foo[["estMTD"]];
        stage1_enrollment = rep(curr_design[["module1"]]$n, n_sim);
        rm(foo);
      } else {
        stop(paste0("Error in entry ",design_list_reorder[k]," of 'design_list': 'module1' is '", curr_design[["module1"]]$name,"' but must be 'crm', '3pl3', 'empiric', or 'fixed'"));
      }

      module1_dat = cbind(array_id = array_id,
                          scenario = dose_outcome_curves[["scenario"]],
                          design = design_list_reorder[k],
                          module1_dat,
                          eff = NA,
                          eff_prob = NA);
      #Generate efficacy data
      module1_dat[,"eff_prob"] = store_eff_curves["stage1",module1_dat[,"dose_num"]];
      module1_dat[,"eff"] = rbinom(nrow(module1_dat),1,module1_dat[,"eff_prob"]);
    }
    store_module1_dat[[k]] = module1_dat;
    store_stage1_estMTD[[k]] = stage1_estMTD;
    store_stage1_enrollment[[k]] = stage1_enrollment;

    #Module 2: First efficacy----

    #Safety-only code dictionary
    #1N = All dose levels found to be unsafe during Stage 1, either during or at end of enrollment
    #1Y = Stage 1 enrollment completed successfully (without regard to efficacy analysis)
    #Safety + Efficacy code dictionary
    #1TN = Trial stopped for safety sometime during stage 1 (so no efficacy analysis was possible)
    #1EN = Trial conducted the stage 1 efficacy analysis and failed
    #1Y = Trial conducted the stage 1 efficacy analysis and passed
    set.seed(seeds_by_module[2]);

    if(all(shared_design_elements[k,c("module1","module2")] != k) &&
       duplicated(shared_design_elements[1:k,paste0("module",1:2)])[k]) {
      if(verbose) {cat("reusing Module 2 conclusions...\n\n");}
      which_duplicate_of =
        which(rowSums(abs(shared_design_elements[1:(k-1),paste0("module",1:2),drop = F] -
                            shared_design_elements[rep(k, k-1),paste0("module",1:2),drop=F])) < eps)[1];
      stage1_summary = store_stage1_summary[[which_duplicate_of]];
      stage1_summary[,"design"] = design_list_reorder[k];
      stage1_RP2D = store_stage1_RP2D[[which_duplicate_of]];
      rm(which_duplicate_of);
    } else {

      n = x = numeric(n_dose);
      names(n) = names(x) = seq_len(n_dose);
      curr_accept_dose = which(store_eff_curves["stage1", grep("_accept",colnames(store_eff_curves))]==1);
      if(near(length(curr_accept_dose), 0)) { curr_accept_dose = 0;}

      stage1_summary = data.frame(cbind(array_id,
                                        dose_outcome_curves[["scenario"]],
                                        design_list_reorder[k],
                                        seq_len(n_sim),
                                        stage1_enrollment,
                                        NA, NA, NA, NA, NA, NA, NA, NA, NA, NA));
      colnames(stage1_summary) = c("array_id",
                                   "scenario",
                                   "design",
                                   "sim_id",
                                   "n_total_enrolled",
                                   "n_possible_enrolled","estMTD","estMTDCode","trueMTD","RP2D","RP2DCode","RP2DAcceptable","bestP2D","estEff_at_estMTD","num_at_estMTD");
      #MTD Summary
      if(curr_design$module1[["name"]] == "3pl3") {
        # We use the convention that the number of possible subjects enrolled
        # in a 3+3 is always defined  to be the *actual* enrollment. This is
        # not perfectly satisfactory but is necessary because for a 3+3 trial
        # to stop due to >1 DLT is different than a CRM stopping because *all*
        # dose levels are  found to be too toxic.
        stage1_summary[,"n_possible_enrolled"] = stage1_enrollment;
      } else {
        stage1_summary[,"n_possible_enrolled"] = curr_design$module1$n;
      }
      stage1_summary[,"estMTD"] = stage1_estMTD;
      stage1_summary[,"estMTDCode"] = ifelse(near(stage1_estMTD, 0), "1N", "1Y");
      stage1_summary[,"trueMTD"] = trueMTD;
      stage1_summary[,"bestP2D"] = max(curr_accept_dose);
      stage1_summary[,"num_at_estMTD"] =
        cbind(module1_dat[,c("sim_id","dose_num")],
              estMTD = stage1_estMTD[module1_dat[,"sim_id"]]) %>%
        mutate(receivedEstMTD = near(dose_num, estMTD)) %>%
        select(sim_id, receivedEstMTD) %>%
        group_by(sim_id) %>%
        summarise(receivedEstMTD = sum(receivedEstMTD)) %>%
        select(receivedEstMTD);

      #Did module1 enroll all patients (and estimate MTD)?
      module1_finished = stage1_estMTD > 0;
      #RP2D will be the MTD if it meets the efficacy target, otherwise it will be 0.

      for(curr_sim in seq_len(n_sim)) {

        if(curr_design[["module2"]]$name == "none") {
          estEff_at_estMTD = NA;
          RP2D = stage1_estMTD[curr_sim];
          RP2DCode = ifelse(near(RP2D, 0), "1TN","1Y");
          RP2DAcceptable = ifelse(near(RP2D, 0),
                                  near(sum(store_eff_curves["stage1",paste0("dose",1:n_dose,"_accept")]), 0),
                                  near(store_eff_curves["stage1",paste0("dose",1:n_dose,"_accept")][RP2D], 1));

        } else if(module1_finished[curr_sim]) {
          # Start fresh
          x = x * 0;
          n = n * 0;

          #RP2D Summary
          curr_data = filter(module1_dat, near(sim_id, curr_sim));
          curr_summary = xtabs(~ dose_num, curr_data);
          n[rownames(curr_summary)] = curr_summary;
          if(any(near(curr_data[,"eff"], 1))) {
            curr_summary = xtabs(~ dose_num, filter(curr_data,near(eff, 1)));
            x[rownames(curr_summary)] = curr_summary;
          }
          rm(curr_data, curr_summary);

          if(curr_design[["module2"]]$name == "bayes_isoreg") {

            #Bayesian isotonic regression via stan
            bayes_iso_fit = bayesian_isotonic(data_grouped =
                                                bind_cols(x = as.numeric(names(x)),
                                                          y = x,
                                                          n = n) %>%
                                                arrange(x),
                                              stan_args = list(
                                                local_dof_stan = 1,
                                                global_dof_stan = 1,
                                                alpha_scale_stan = curr_design[["module2"]]$alpha_scale,
                                                slab_precision_stan = 1),
                                              conf_level = curr_design[["module2"]]$prob_threshold,
                                              conf_level_direction = "lower",
                                              n_mc_warmup = stan_args$n_mc_warmup,
                                              n_mc_samps = stan_args$n_mc_samps,
                                              mc_chains = stan_args$mc_chains,
                                              mc_thin = stan_args$mc_thin,
                                              mc_stepsize = stan_args$mc_stepsize,
                                              mc_adapt_delta = stan_args$mc_adapt_delta,
                                              mc_max_treedepth = stan_args$mc_max_treedepth,
                                              ntries = stan_args$ntries);

            estEff_at_estMTD =
              bayes_iso_fit$data_grouped %>%
              filter(x == stage1_estMTD[curr_sim]) %>%
              pull(model_mean_prob)
            if(bayes_iso_fit$data_grouped %>%
               filter(x == stage1_estMTD[curr_sim]) %>%
               pull(model_lower_ci_prob) >=
               primary_objectives["eff_target"]) {
              RP2D = stage1_estMTD[curr_sim];
            } else {
              RP2D = 0;
            }
            rm(bayes_iso_fit);

            #Independent beta priors
          } else if(curr_design[["module2"]]$name == "bayes") {
            beta_shape1 = x[stage1_estMTD[curr_sim]] + curr_design[["module2"]]$prior_n_per * curr_design[["module2"]]$prior_mean[stage1_estMTD[curr_sim]];
            beta_shape2 = (n[stage1_estMTD[curr_sim]] - x[stage1_estMTD[curr_sim]]) + curr_design[["module2"]]$prior_n_per * (1 - curr_design[["module2"]]$prior_mean[stage1_estMTD[curr_sim]]);
            estEff_at_estMTD = beta_shape1 / (beta_shape1 + beta_shape2);
            RP2D = ifelse(qbeta(curr_design[["module2"]]$prob_threshold, beta_shape1, beta_shape2,lower.tail=F) >= primary_objectives["eff_target"], stage1_estMTD[curr_sim], 0);
            rm(beta_shape1,beta_shape2);
          } else if(curr_design[["module2"]]$name == "inverted_score") {
            #Assume a response rate of 0 / 0.5 if no patients have actually been enrolled to the current estimated MTD
            if(n[stage1_estMTD[curr_sim]] > 0) {
              foo = binom.wilson(x = x[stage1_estMTD[curr_sim]],n = n[stage1_estMTD[curr_sim]],conf.level = 1 - 2*(1-curr_design[["module2"]]$ci_level_onesided));
              estEff_at_estMTD = foo[,"mean"];
            } else {
              foo = binom.wilson(x = 0, n = 0.5, conf.level = 1 - 2*(1-curr_design[["module2"]]$ci_level_onesided));
              estEff_at_estMTD = NA;
            }
            RP2D = ifelse(foo[,"upper"] >= primary_objectives["eff_target"], stage1_estMTD[curr_sim], 0);
            rm(foo);
          } else if(curr_design[["module2"]]$name == "min_num_resp") {
            estEff_at_estMTD = NA;
            RP2D = ifelse(x[stage1_estMTD[curr_sim]] >= curr_design[["module2"]]$number, stage1_estMTD[curr_sim], 0);
          } else if(curr_design[["module2"]]$name == "min_pct_resp") {
            estEff_at_estMTD = NA;
            RP2D = ifelse((n[stage1_estMTD[curr_sim]] > 0) && (x[stage1_estMTD[curr_sim]]/n[stage1_estMTD[curr_sim]] >= curr_design[["module2"]]$percent), stage1_estMTD[curr_sim], 0);
          } else {
            stop(paste0("Error in entry ",design_list_reorder[k]," of 'design_list': 'module2' is '", curr_design[["module2"]]$name,"' but  must be 'none', 'bayes', 'bayes_isoreg', or 'inverted_score'"));
          }

          RP2DCode = ifelse(near(RP2D, 0), "1EN","1Y");
          RP2DAcceptable = ifelse(near(RP2D, 0),
                                  sum(store_eff_curves["stage1",paste0("dose",1:n_dose,"_accept")]) == 0,
                                  store_eff_curves["stage1",paste0("dose",1:n_dose,"_accept")][RP2D] == 1);

        } else {
          estEff_at_estMTD = NA;
          RP2D = 0;
          RP2DCode = "1TN"
          RP2DAcceptable = sum(store_eff_curves["stage1",paste0("dose",1:n_dose,"_accept")]) == 0;

        }
        stage1_summary[curr_sim,"RP2D"] = RP2D;
        stage1_summary[curr_sim,"RP2DCode"] = RP2DCode;
        stage1_summary[curr_sim,"RP2DAcceptable"] = RP2DAcceptable;
        stage1_summary[curr_sim,"estEff_at_estMTD"] = estEff_at_estMTD;
        rm(RP2D,RP2DCode,RP2DAcceptable,estEff_at_estMTD);
      }

      stage1_RP2D = as.numeric(stage1_summary[,"RP2D"]);
      store_stage1_summary[[k]] = stage1_summary;
      store_stage1_RP2D[[k]] = stage1_RP2D;
      rm(curr_accept_dose,curr_sim,module1_finished,n,x);
    }

    #Module 3: Second Dose Assignment----
    set.seed(seeds_by_module[3]);
    #This is the case the design is completely identical to a previously simulated design up to this point
    if(curr_design[["module3"]]$name != "none" &&
       all(shared_design_elements[k,c("module1","module2","module3")] != k) &&
       duplicated(shared_design_elements[1:k,paste0("module",1:3)])[k]) {
      if(verbose) {cat("reusing Module 3 conclusions...\n\n");}
      which_duplicate_of =
        which(rowSums(abs(shared_design_elements[1:(k-1),paste0("module",1:3),drop = F] -
                            shared_design_elements[rep(k, k-1),paste0("module",1:3),drop = F])) < eps)[1];
      module3_dat = store_module3_dat_postmodule2[[which_duplicate_of]];
      if(nrow(module3_dat) > 0) {
        module3_dat[,"design"] = design_list_reorder[k];
      }
      stage2_estMTD = store_stage2_estMTD_postmodule2[[which_duplicate_of]];
      stage2_enrollment = store_stage2_enrollment_postmodule2[[which_duplicate_of]];
      rm(which_duplicate_of);
    } else {
      if(curr_design[["module3"]]$name == "continue_crm") {
        module3_dat =
          module3_dat  %>%
          filter(sim_id %in% which(stage1_RP2D >0));
        stage2_estMTD = (stage1_RP2D > 0) * stage2_estMTD;
        stage2_enrollment = (stage1_RP2D > 0) * stage2_enrollment;
      } else if(any(stage1_RP2D > 0)) {
        if(curr_design[["module3"]]$name == "crm") {

          titecrm_args = starting_titecrm_args;
          titecrm_args$n = curr_design[["module3"]]$n;
          titecrm_args$cohort.size = curr_design[["module3"]]$dose_cohort_size;
          titecrm_args$first.cohort.only = curr_design[["module3"]]$dose_cohort_size_first_only;
          titecrm_args$earliest_stop = curr_design[["module3"]]$earliest_stop;

          foo = sim_twostage_crm(n_sim = n_sim,
                                 titecrm_args = titecrm_args,
                                 second_stage_start_after = Inf,
                                 first_stage_label = 2,
                                 sim_specific_start_id = stage1_enrollment + 1,
                                 sim_specific_prior = matrix(rep(curr_design[["module3"]]$skeleton, n_sim), nrow=n_sim,
                                                             byrow=TRUE),
                                 sim_specific_x0 = stage1_RP2D,
                                 sim_specific_scale = rep(curr_design[["module3"]]$beta_scale, n_sim));

          module3_dat = data.frame(foo[["all_results"]]);
          stage2_estMTD = foo[["first_stage_estMTD"]];
          stage2_enrollment = foo[["first_stage_enrollment"]];
          rm(foo);

        } else if(curr_design[["module3"]]$name == "3pl3") {
          foo = sim_3pl3(n_sim = n_sim,
                         true_tox_curve = dose_outcome_curves[["tox_curve"]],
                         stage_label = 2,
                         sim_specific_start_id = stage1_enrollment + 1,
                         sim_specific_dose_start = stage1_RP2D);

          module3_dat = data.frame(foo[["all_results"]]);
          stage2_estMTD = foo[["estMTD"]];
          stage2_enrollment = foo[["enrollment"]];
          rm(foo);
        } else if(curr_design[["module3"]]$name == "empiric") {

          foo = sim_empiric_dec(n_sim = n_sim,
                                true_tox_curve = dose_outcome_curves[["tox_curve"]],
                                stage_label = 2,
                                sim_specific_start_id = stage1_enrollment + 1,
                                sim_specific_dose_start = stage1_RP2D,
                                max_n_per_dec = curr_design[["module3"]]$n,
                                module_rule = curr_design[["module3"]]$rule,
                                thresh_decrease = curr_design[["module3"]]$thresh_decrease,
                                first_patient_look = curr_design[["module3"]]$first_patient_look);
          stage2_estMTD = foo[["estMTD"]];
          stage2_enrollment = foo[["enrollment"]];
          module3_dat = data.frame(foo[["all_results"]]);
          rm(foo);
        } else if(curr_design[["module3"]]$name == "fixed") {
          foo = sim_empiric_dec(n_sim = n_sim,
                                true_tox_curve = dose_outcome_curves[["tox_curve"]],
                                stage_label = 2,
                                sim_specific_start_id = stage1_enrollment + 1,
                                sim_specific_dose_start = stage1_RP2D,
                                max_n_per_dec = curr_design[["module3"]]$n,
                                module_rule = "local",
                                thresh_decrease = Inf);
          stage2_estMTD = foo[["estMTD"]];
          stage2_enrollment = foo[["enrollment"]];
          module3_dat = data.frame(foo[["all_results"]]);
          rm(foo);
        } else if(curr_design[["module3"]]$name == "none") {
          module3_dat = data.frame(matrix(NA,nrow = 0, ncol = 11));
          colnames(module3_dat) = c("array_id",
                                    "scenario",
                                    "design",
                                    "sim_id",
                                    "subj_id",
                                    "dose_num",
                                    "tox_prob",
                                    "tox",
                                    "stage",
                                    "eff",
                                    "eff_prob");
          stage2_estMTD = stage1_estMTD;
          stage2_enrollment = numeric(n_sim);
        } else {
          stop(paste0("Error in entry ",design_list_reorder[k]," of 'design_list': 'module3' is '", curr_design[["module1"]]$name,"' but must be 'crm', 'continue_crm', '3pl3', 'empiric', 'fixed', or 'none'"));
        }

        if(any(near(stage1_RP2D, 0))) {
          stage2_estMTD[which(near(stage1_RP2D, 0))] = 0;
        }

        if(nrow(module3_dat) > 0) {
          module3_dat = cbind(array_id = array_id,
                              scenario = dose_outcome_curves[["scenario"]],
                              design = design_list_reorder[k],
                              module3_dat,
                              eff = NA,
                              eff_prob = NA);

          #Generate efficacy data
          module3_dat[,"eff_prob"] = store_eff_curves["stage1",module3_dat[,"dose_num"]];
          module3_dat[,"eff"] = rbinom(nrow(module3_dat),1,module3_dat[,"eff_prob"]);
        }
      } else if(curr_design[["module3"]]$name %in% c("continue_crm","crm","3pl3","empiric","fixed","none")) {
        #The module is valid but there is no need to simulate data because none of the simulated trials made it to stage 2
        module3_dat = data.frame(matrix(NA,nrow = 0, ncol = 11));
        colnames(module3_dat) = c("array_id",
                                  "scenario",
                                  "design",
                                  "sim_id",
                                  "subj_id",
                                  "dose_num",
                                  "tox_prob",
                                  "tox",
                                  "stage",
                                  "eff",
                                  "eff_prob");
        stage2_estMTD = numeric(n_sim); #Is there some reason we're setting this to n_sim? 0 seems more intuitive to me.
        stage2_enrollment = numeric(n_sim);
      } else {
        stop(paste0("Error in entry ",design_list_reorder[k]," of 'design_list': 'module3' is '", curr_design[["module1"]]$name,"' but must be 'crm', 'continue_crm', '3pl3', 'empiric', 'fixed', or 'none'"));
      }

      store_module3_dat_postmodule2[[k]] = module3_dat;
      store_stage2_estMTD_postmodule2[[k]] = stage2_estMTD;
      store_stage2_enrollment_postmodule2[[k]] = stage2_enrollment;

    }
    #Module 5: Second efficacy----

    #Safety-only code dictionary
    #1TN = Stage 2 did not open (all dose levels found to be unsafe during Stage 1, either during or at end of enrollment)
    #1EN = Stage 2 did not open (trial conducted the stage 1 efficacy analysis and failed)
    #2TN = Stage 2 opened but all dose levels were found to be unsafe during Stage 2, either during or at end of enrollment
    #2Y = Stage 2 opened and successfully estimated an MTD (without regard to efficacy at end of trial)
    #Safety + Efficacy code dictionary
    #1TN = Stage 2 did not open (trial stopped for safety sometime during stage 1, either during or at end of enrollment)
    #1EN = Stage 2 did not open (trial conducted the stage 1 efficacy analysis and failed)
    #2TN = Stage 2 opened but all dose levels were found to be unsafe during Stage 2, either during or at end of enrollment
    #2EN = Trial conducted the stage 2 efficacy analysis and failed
    #2Y = Trial conducted the stage 2 efficacy analysis and passed
    set.seed(seeds_by_module[4]);
    if(curr_design[["module3"]]$name == "none") {
      curr_design[["module3"]]$n = 0;
    }
    n = x = numeric(n_dose);
    names(n) = names(x) = seq_len(n_dose);
    curr_accept_dose = which(near(store_eff_curves["stage2", grep("_accept",colnames(store_eff_curves))], 1));
    if(near(length(curr_accept_dose), 0)) { curr_accept_dose = 0;}

    stage2_summary = data.frame(cbind(array_id,
                                      dose_outcome_curves[["scenario"]],
                                      design_list_reorder[k],
                                      seq_len(n_sim),
                                      stage1_enrollment + stage2_enrollment,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA));
    colnames(stage2_summary) = c("array_id",
                                 "scenario",
                                 "design",
                                 "sim_id",
                                 "n_total_enrolled",
                                 "n_possible_enrolled","estMTD","estMTDCode","trueMTD","RP2D","RP2DCode","RP2DAcceptable","bestP2D","estEff_at_estMTD","num_at_estMTD");
    #MTD Summary
    if(curr_design$module3[["name"]] == "3pl3") {
      stage2_summary[,"n_possible_enrolled"] = stage1_summary[,"n_possible_enrolled"] + stage2_enrollment;
    } else {
      stage2_summary[,"n_possible_enrolled"] = stage1_summary[,"n_possible_enrolled"] + curr_design[["module3"]]$n;
    }
    stage2_summary[,"estMTD"] = stage2_estMTD;
    stage2_summary[,"estMTDCode"] =
      ifelse(stage1_summary[,"estMTDCode"] == "1N", "1TN",
             ifelse(stage1_summary[,"RP2DCode"] == "1EN", "1EN",
                    ifelse(near(stage2_estMTD, 0), "2TN", "2Y")));
    stage2_summary[,"trueMTD"] = trueMTD;
    stage2_summary[,"bestP2D"] = max(curr_accept_dose);
    stage2_summary[,"num_at_estMTD"] =
      rbind(cbind(module1_dat[,c("sim_id","dose_num")], estMTD = stage2_estMTD[module1_dat[,"sim_id"]]),
            cbind(module3_dat[,c("sim_id","dose_num")], estMTD = stage2_estMTD[module3_dat[,"sim_id"]])) %>%
      arrange(sim_id) %>%
      mutate(receivedEstMTD = near(dose_num, estMTD)) %>%
      select(sim_id, receivedEstMTD) %>%
      group_by(sim_id) %>%
      summarise(receivedEstMTD = sum(receivedEstMTD)) %>%
      select(receivedEstMTD);

    #Did module 3 start (!is.na) and finish?
    module3_finished = (stage2_estMTD > 0);
    #RP2D will be the MTD if it meets the efficacy target, otherwise it will be 0.

    for(curr_sim in seq_len(n_sim)) {

      if(curr_design[["module4"]]$name == "none") {
        #If both module3 and module4 were not run, then trial really ended at module 2 and those results should be carried forward
        if(curr_design[["module3"]]$name == "none") {
          estEff_at_estMTD = stage1_summary[curr_sim,"estEff_at_estMTD"];
        } else {
          estEff_at_estMTD = NA;
        }
        # These are left the same regardless. The 'RP2DCode' is an issue of semantics, when did the recommendation occur:
        # at the second stage that didn't happen or at the first stage?
        RP2D = stage2_estMTD[curr_sim];
        RP2DCode = stage2_summary[curr_sim,"estMTDCode"];
        RP2DAcceptable = ifelse(near(RP2D, 0),
                                near(sum(store_eff_curves["stage2",paste0("dose",1:n_dose,"_accept")]), 0),
                                near(store_eff_curves["stage2",paste0("dose",1:n_dose,"_accept")][RP2D], 1));

      } else if(module3_finished[curr_sim]) {
        # Start fresh
        x = x * 0;
        n = n * 0;

        if(curr_design[["module4"]][["include_stage1_data"]]) {
          curr_dat = rbind(filter(module1_dat, near(sim_id, curr_sim)),
                           filter(module3_dat, near(sim_id, curr_sim)));
        } else {
          curr_dat = filter(module3_dat, near(sim_id, curr_sim));
        }

        #RP2D Summary
        curr_summary = xtabs(~ dose_num, curr_dat);
        n[rownames(curr_summary)] = curr_summary;
        if(any(near(curr_dat[,"eff"], 1))) {
          curr_summary = xtabs(~ dose_num, filter(curr_dat, near(eff, 1)));
          x[rownames(curr_summary)] = curr_summary;
        }
        rm(curr_summary,curr_dat);

        if(curr_design[["module4"]]$name == "bayes_isoreg") {

          #Bayesian isotonic regression via stan
          bayes_iso_fit = bayesian_isotonic(data_grouped =
                                              bind_cols(x = as.numeric(names(x)),
                                                        y = x,
                                                        n = n) %>%
                                              arrange(x),
                                            stan_args = list(
                                              local_dof_stan = 1,
                                              global_dof_stan = 1,
                                              alpha_scale_stan = curr_design[["module4"]]$alpha_scale,
                                              slab_precision_stan = 1),
                                            conf_level = curr_design[["module4"]]$prob_threshold,
                                            conf_level_direction = "lower",
                                            n_mc_warmup = stan_args$n_mc_warmup,
                                            n_mc_samps = stan_args$n_mc_samps,
                                            mc_chains = stan_args$mc_chains,
                                            mc_thin = stan_args$mc_thin,
                                            mc_stepsize = stan_args$mc_stepsize,
                                            mc_adapt_delta = stan_args$mc_adapt_delta,
                                            mc_max_treedepth = stan_args$mc_max_treedepth,
                                            ntries = stan_args$ntries);

          #estEff_at_estMTD = bayes_iso_fit$posterior_means[stage2_estMTD[curr_sim]];
          #RP2D =
          #  ifelse(bayes_iso_fit$posterior_intervals["lower",stage2_estMTD[curr_sim]] >= primary_objectives["eff_target"], stage2_estMTD[curr_sim],0);
          estEff_at_estMTD =
            bayes_iso_fit$data_grouped %>%
            filter(x == stage2_estMTD[curr_sim]) %>%
            pull(model_mean_prob)
          if(bayes_iso_fit$data_grouped %>%
             filter(x == stage2_estMTD[curr_sim]) %>%
             pull(model_lower_ci_prob) >=
             primary_objectives["eff_target"]) {
            RP2D = stage2_estMTD[curr_sim];
          } else {
            RP2D = 0;
          }
          rm(bayes_iso_fit);
          #Independent beta priors
        } else if(curr_design[["module4"]]$name == "bayes") {
          beta_shape1 = x[stage2_estMTD[curr_sim]] + curr_design[["module4"]]$prior_n_per * curr_design[["module4"]]$prior_mean[stage2_estMTD[curr_sim]];
          beta_shape2 = (n[stage2_estMTD[curr_sim]] - x[stage2_estMTD[curr_sim]]) + curr_design[["module4"]]$prior_n_per * (1 - curr_design[["module4"]]$prior_mean[stage2_estMTD[curr_sim]]);
          estEff_at_estMTD = beta_shape1 / (beta_shape1 + beta_shape2);
          RP2D = ifelse(qbeta(curr_design[["module4"]]$prob_threshold, beta_shape1, beta_shape2,lower.tail = F) >= primary_objectives["eff_target"], stage2_estMTD[curr_sim], 0);
        } else if(curr_design[["module4"]]$name == "inverted_score") {
          #1-sided confidence interval
          if(n[stage2_estMTD[curr_sim]] > 0) {
            foo = binom.wilson(x = x[stage2_estMTD[curr_sim]],n = n[stage2_estMTD[curr_sim]],conf.level = 1 - 2*(1-curr_design[["module4"]]$ci_level_onesided));
            estEff_at_estMTD = foo[,"mean"];
            RP2D = ifelse(foo[,"lower"] >= primary_objectives["eff_target"], stage2_estMTD[curr_sim], 0);
            rm(foo);
          } else {
            estEff_at_estMTD = NA;
            RP2D = 0;
          }
        } else if(curr_design[["module4"]]$name == "min_num_resp") {
          estEff_at_estMTD = NA;
          RP2D = ifelse(x[stage2_estMTD[curr_sim]] >= curr_design[["module4"]]$number, stage2_estMTD[curr_sim], 0);
        } else if(curr_design[["module4"]]$name == "min_pct_resp") {
          estEff_at_estMTD = NA;
          RP2D = ifelse((n[stage2_estMTD[curr_sim]] > 0) && (x[stage2_estMTD[curr_sim]]/n[stage2_estMTD[curr_sim]] >= curr_design[["module4"]]$percent), stage2_estMTD[curr_sim], 0);
        } else {
          stop(paste0("Error in entry ",design_list_reorder[k]," of 'design_list': 'module4' is '", curr_design[["module4"]]$name,"' but  must be 'none', 'bayes', 'bayes_isoreg', or 'inverted_score'"));
        }

        RP2DCode = ifelse(near(RP2D, 0), "2EN","2Y");
        RP2DAcceptable = ifelse(near(RP2D, 0),
                                near(sum(store_eff_curves["stage2",paste0("dose",1:n_dose,"_accept")]), 0),
                                near(store_eff_curves["stage2",paste0("dose",1:n_dose,"_accept")][RP2D], 1));

      } else {
        estEff_at_estMTD = NA;
        RP2D = 0;
        RP2DCode = stage2_summary[curr_sim,"estMTDCode"];
        RP2DAcceptable = near(sum(store_eff_curves["stage2",paste0("dose",1:n_dose,"_accept")]), 0);
      }
      stage2_summary[curr_sim,"RP2D"] = RP2D;
      stage2_summary[curr_sim,"RP2DCode"] = RP2DCode;
      stage2_summary[curr_sim,"RP2DAcceptable"] = RP2DAcceptable;
      stage2_summary[curr_sim,"estEff_at_estMTD"] = estEff_at_estMTD;

      rm(RP2D,RP2DCode,RP2DAcceptable,estEff_at_estMTD);
    }
    stage2_RP2D = stage2_summary[,"RP2D"];
    rm(curr_accept_dose,curr_sim,module3_finished,n,x);

    #End of modules----

    sim_data_stage1 = rbind(sim_data_stage1, stage1_summary);
    sim_data_stage2 = rbind(sim_data_stage2, stage2_summary);

    #The next line joins both module1_dat and module3_dat with stage2_summary.
    #That module1_dat is joined with stage2_summary (and not stage1_summary) is not a typo;
    #rather, for consistency, the patient level data is always joined with
    #the outcome at the very end of the trial (or when the trial stopped).
    curr_patient_data =
      arrange(rbind(
        left_join(module1_dat,
                  select(stage2_summary,-array_id,-scenario,-design), by = "sim_id"),
        left_join(module3_dat,
                  select(stage2_summary,-array_id,-scenario,-design), by = "sim_id")),
        sim_id,subj_id);
    #Fill in any zero dose assignments due to the trial having stopped prior to its maximum possible enrollment
    if(any(curr_patient_data[,"n_possible_enrolled"] - curr_patient_data[,"n_total_enrolled"] > 0)) {
      which_sims_to_pad = unique(curr_patient_data[which(curr_patient_data[,"n_possible_enrolled"] - curr_patient_data[,"n_total_enrolled"] > 0),"sim_id"]);
      for(curr_sim in which_sims_to_pad) {
        template_row = filter(curr_patient_data, near(sim_id, curr_sim), near(subj_id, 1));
        template_row[,c("dose_num","tox_prob","tox","eff","eff_prob")] = 0;
        template_row = template_row[rep(1,template_row[["n_possible_enrolled"]] - template_row[["n_total_enrolled"]]),];
        template_row[,"subj_id"] = (1:nrow(template_row)) + filter(curr_patient_data, near(sim_id, curr_sim)) %>% select(subj_id) %>% max();
        template_row[,"stage"] = 1 + (template_row[,"subj_id"] > filter(stage1_summary, near(sim_id, curr_sim))[["n_possible_enrolled"]]);
        curr_patient_data = rbind(curr_patient_data, template_row);
        rm(template_row);
      }
      rm(which_sims_to_pad, curr_sim);
      curr_patient_data =
        curr_patient_data %>%
        arrange(sim_id,subj_id)
    }

    patient_data_to_return =
      rbind(patient_data_to_return, curr_patient_data);

    rm(module1_dat, module3_dat,curr_patient_data,
       stage1_enrollment,stage1_estMTD,stage1_RP2D,stage1_summary,
       stage2_enrollment,stage2_estMTD,stage2_RP2D,stage2_summary);
  }

  patient_data_to_return =
    patient_data_to_return %>%
    arrange(design, sim_id, subj_id) %>%
    mutate(sim_id = sim_labels[sim_id],
           design = design_labels[design]);

  sim_data_stage1 =
    sim_data_stage1 %>%
    arrange(design, sim_id) %>%
    mutate(sim_id = sim_labels[sim_id],
           design = design_labels[design]);

  sim_data_stage2 =
    sim_data_stage2 %>%
    arrange(design, sim_id) %>%
    mutate(sim_id = sim_labels[sim_id],
           design = design_labels[design]);

  design_list = design_list[order(design_list_reorder)];
  design_description = design_description[order(design_list_reorder),];
  shared_design_elements = shared_design_elements[order(design_list_reorder),];

  list(patient_data = patient_data_to_return,
       sim_data_stage1 = sim_data_stage1,
       sim_data_stage2 = sim_data_stage2,
       dose_outcome_curves = dose_outcome_curves,
       titecrm_args = starting_titecrm_args,
       design_list = design_list,
       design_description = design_description,
       shared_design_elements = shared_design_elements,
       random_seed = random_seed);
}

