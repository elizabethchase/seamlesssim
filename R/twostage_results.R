#' This function processes the output from twostage_simulator to produce several
#' basic plots and summary tables.
#'
#' For more information on function inputs and features, please see the vignette
#' and Boonstra (2020) on Arxiv.
#'
#' @param csv A logical value indicating whether the output from twostage_simulator has been saved
#' as .csv file(s) (if .csv, then TRUE).
#' @param stage2folder A character string giving the filepath to a folder containing
#' nothing but sim_data_stage2 output from twostage_simulator, saved as .csv file(s).
#' This should only be provided if csv = TRUE.
#' @param patientdatfolder A character string giving the filepath to a folder containing
#' nothing but patient_data output from twostage_simulator, saved as .csv file(s).
#' This should only be provided if csv = TRUE.
#' @param dose_outcome_curves_list A list of lists, with each list element containing three named
#' elements and an optional fourth
#' element: tox_curve, eff_curve, scenario, and, optionally, eff_curve_stage2. tox_curve is
#' the true toxicity curve of the doses; eff_curve is the true efficacy curve of the doses;
#' scenario is an identifier of which true data-generating scenario is being run (meant to be
#' helpful to the user when calling this function multiple times for different data-generating
#' scenarios). This should be the same list that was used in twostage_simulator to generate all
#' the simulated trials. This should only be provided if csv = TRUE.
#' @param files A character string giving the filepath to a folder containing
#' nothing but the raw output from twostage_simulator, saved as .Rds file(s). Alternatively, this can
#' a character vector containing the names of all of the raw output from twostage_simulator, stored
#' in the workspace. This should only be provided if csv = FALSE.
#' @param filepath A logical value indicating whether files contains a filepath (filepath = TRUE) or
#' a character vector of file names stored in the workspace (filepath = FALSE). The default is TRUE.
#' This should only be provided if csv = FALSE.
#' @param primary_objectives A list containing three named elements: tox_target, tox_delta_no_exceed,
#' and eff_target, such that tox_target is between 0 and 1, tox_delta_no_exceed is between 0 and
#' (1 - tox_target), and eff_target is between 0 and 1. This should be the same list that was used in
#' twostage_simulator to generate the simulated trials.
#' @param design_labels A character vector giving labels for the different designs. If left NULL,
#' the number of the design will be used as its label. This should be the same length and the same order
#' as design_list inputted in twostage_simulator.
#' @param scen_per_page A numeric value indicating the number of data-generating scenarios that should
#' be printed per page. This can be at most 10 (the default is 10).
#' @param design_per_page A numeric value indicating the number of designs that should be printed per
#' page. This can be at most 3 (the default is 3).
#' @return The function returns a named list containing two items: "plots" and "tables".
#' tables contains:
#' \describe{
#'    \item{generating_params_for_display}{This is a table giving the true efficacy and toxicity of
#'    each dose for each scenario considered.}
#'    \item{acc_dose_rec_table}{This is a table giving which proportion of trials for each design and
#'    scenario combination recommended an acceptable dose (a dose meeting the toxicity and efficacy
#'    standards, if not the best dose) at the end of the trial.}
#'    \item{mean_patients_table}{This is a table giving the mean number of patients enrolled for each
#'    design and scenario combination across the different simulated trials.}
#' }
#' plots contains:
#' \describe{
#'    \item{gen_param_plot}{This is a list of plots providing the same information as in generating_params_for_display.
#'    If the number of designs under consideration exceeds design_per_page or the number of
#'    scenarios under consideration exceeds scen_per_page, the list will contain multiple plots. To call the
#'    first plot, we would run: results$plots$gen_param_plot[[1]] and so on for further plots.}
#'    \item{acc_dose_rec_plot}{This is a plot giving which proportion of trials for each design and
#'    scenario combination recommended an acceptable dose, unacceptable dose, or made no recommendation
#'    at all. If the number of designs under consideration exceeds design_per_page or the number of
#'    scenarios under consideration exceeds scen_per_page, the list will contain multiple plots. To call the
#'    first plot, we would run: results$plots$acc_dose_rec_plot[[1]] and so on for further plots.}
#'    \item{dose_over_time_plot}{This is a filled barplot giving the distribution of dose assignments by
#'    time, where time is measured by patient number. To call the
#'    first plot, we would run: results$plots$dose_over_time_plot[[1]] and so on for further plots.}
#'    \item{n_patients_plot}{This is a boxplot of the number of patients enrolled for each design and
#'    scenario combination across the different simulated trials.}
#'    \item{n_patients_RP2D_plot}{This is a boxplot of the number of patients who received the final
#'    recommended dose for each design and scenario combination across the different simulated trials.
#'    Note that the final recommended dose may not be safe or effective--this is merely a measure of
#'    how much patient data the design will yield for the final dose it recommends.}
#'    \item{prop_patients_accep_plot}{This is a boxplot of the number of patients who received an
#'    acceptable dose (a dose meeting the toxicity and efficacy standards, if not the best dose) for
#'    each design and scenario combination across the different simulated trials.}
#' }
#' @importFrom sjmisc is_empty
#' @importFrom dplyr %>% mutate arrange near group_by tally n ungroup select filter
#' summarize left_join
#' @importFrom tidyr gather spread
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_grid labs scale_y_continuous expansion
#' scale_color_manual guides guide_legend theme element_text margin element_blank geom_bar geom_text
#' geom_col scale_y_reverse scale_x_continuous scale_fill_manual scale_size_manual geom_boxplot
#' position_dodge2 ylab xlab scale_fill_discrete unit
#' @importFrom forcats fct_inorder fct_relabel
#' @import RColorBrewer
#' @export
twostage_results <- function(csv = FALSE,
                             stage2folder = NULL,
                             patientdatfolder = NULL,
                             dose_outcome_curves_list = NULL,
                             files = NULL,
                             filepath = TRUE,
                             primary_objectives = NULL,
                             design_labels = NULL,
                             scen_per_page = 10,
                             design_per_page = 3){

  if (scen_per_page > 10){
    stop("'scen_per_page' can be at most 10")
  }
  if (design_per_page > 3){
    stop("'design_per_page' can be at most 3")
  }
  if (!("tox_target" %in% names(primary_objectives))){
    stop("'primary_objectives' must be a vector with element 'tox_target'")
  }
  if (!("tox_delta_no_exceed" %in% names(primary_objectives))){
    stop("'primary_objectives' must be a vector with element 'tox_delta_no_exceed'")
  }
  if (!("eff_target" %in% names(primary_objectives))){
    stop("'primary_objectives' must be a vector with element 'eff_target'")
  }

  scenario <- NULL
  dose_num <- NULL
  true_dlt_prob <- NULL
  true_eff_prob <- NULL
  type <- NULL
  is_acceptable <- NULL
  design <- NULL
  array_id <- NULL
  sim_id <- NULL
  RP2DCode <- NULL
  RP2DAcceptable <- NULL
  key <- NULL
  value <- NULL
  RP2DCode_truth <- NULL
  set_designation <- NULL
  design_label <- NULL
  bestP2D <- NULL
  sum_n <- NULL
  prop_n <- NULL
  prob <- NULL
  is_acceptable_by_dose_num <- NULL
  text_height <- NULL
  n_total_enrolled <- NULL
  mean_patients <- NULL
  ataccept <- NULL
  total <- NULL
  atRP2D <- NULL
  prop_acc <- NULL

  text_size = 4;
  legend_text_size = 7;
  min_prop_to_write = 0.25;

  generating_params_for_display =
    #matrix(nrow = 0, ncol = 5, dimnames = list(NULL, c("Scenario","True DLT Probability","True Efficacy Probability","True MTD","Acceptable/Desirable Dose")));
    generating_params =
    tibble();

  scen <- c()

  trial_summary = NULL;
  patient_summary <- NULL

  if (!csv){
    # Indicate that results will be read in from .Rds files

    if (typeof(files) != "character") {
      stop("'files' must be a character string or vector");
    }
    if (filepath){
      myfiles <- list.files(path = files, pattern = ".Rds", full.names=TRUE)
      if (is_empty(myfiles)){
        stop("'files' contains no files with extension .Rds OR 'files' does not point to a valid path")
      }
    } else{
      myfiles <- files
    }
    if(!is.null(dose_outcome_curves_list)) {
      warning("'dose_outcome_curves_list' was provided but was ignored; the values of
              each dose-outcome curve were taken from the .Rds files.")
    }

    for (i in 1:length(myfiles)) {
      if (filepath){
        dat <- readRDS(myfiles[i])
      } else{
        dat <- get(myfiles[i])
      }
      trial_summary <- bind_rows(trial_summary, dat$sim_data_stage2)
      patient_summary <- bind_rows(patient_summary, dat$patient_data)
      if (!dat$dose_outcome_curves$scenario %in% scen){
        newrow =
          tibble(Scenario = dat$dose_outcome_curves$scenario,
                 `True DLT Probability` =  paste0("{",paste0(dat$dose_outcome_curves$tox_curve,collapse=","),"}"),
                 `True Efficacy Probability` = paste0("{",paste0(dat$dose_outcome_curves$eff_curve,collapse=","),"}"),
                 `True MTD` =
                   max(c(which(dat$dose_outcome_curves$tox_curve <= primary_objectives[["tox_target"]] +
                                 primary_objectives[["tox_delta_no_exceed"]]),0)));

        curr_admiss = (dat$dose_outcome_curves$eff_curve >= primary_objectives[["eff_target"]]) &
          (dat$dose_outcome_curves$tox_curve <=  primary_objectives[["tox_target"]] +
             primary_objectives[["tox_delta_no_exceed"]]);
        if(sum(curr_admiss) == 0) {
          newrow <-
            mutate(newrow, `Acceptable/Desirable Dose` = as.character(0))
          curr_admiss = c(T, curr_admiss);
        } else if(sum(curr_admiss) == 1) {
          newrow <-
            mutate(newrow, `Acceptable/Desirable Dose` = as.character(which(curr_admiss)))
          curr_admiss = c(F, curr_admiss);
        } else {
          newrow <-
            mutate(newrow, `Acceptable/Desirable Dose` = paste0("{",paste0(which(curr_admiss),collapse=","),"}"))
          curr_admiss = c(F, curr_admiss);
        }

        generating_params_for_display <-
          bind_rows(generating_params_for_display, newrow);

        mtd_as_logical = (0:length(dat$dose_outcome_curves$tox_curve)) == max(c(which(dat$dose_outcome_curves$tox_curve <= primary_objectives[["tox_target"]] +
                                                                                         primary_objectives[["tox_delta_no_exceed"]]),0))
        generating_params =
          bind_rows(generating_params,
                    tibble(
                      scenario = dat$dose_outcome_curves$scenario,
                      dose_num = 0:length(dat$dose_outcome_curves$tox_curve),
                      true_dlt_prob = c(0, dat$dose_outcome_curves$tox_curve),
                      true_eff_prob = c(0, dat$dose_outcome_curves$eff_curve),
                      tox_target = primary_objectives["tox_target"] + primary_objectives["tox_delta_no_exceed"],
                      eff_target = primary_objectives["eff_target"],
                      is_mtd = mtd_as_logical,
                      is_acceptable = curr_admiss))

        scen <- c(scen, dat$dose_outcome_curves$scenario)
      }
      rm(list=c("dat"))
      cat(i,"\n");
    }
  } else{
    if (typeof(stage2folder) != "character") {
      stop("'stage2folder' must be a character string");
    }
    if (typeof(patientdatfolder) != "character") {
      stop("'patientdatfolder' must be a character string");
    }
    stage2files <- list.files(path = stage2folder, pattern = "Stage2Data\\d*.csv", full.names=TRUE)
    if (is_empty(stage2files)){
      stop("'stage2folder' contains no files with extension .csv OR 'stage2folder' does not point to a valid path")
    }

    patientdatfiles <- list.files(path = patientdatfolder, pattern = "PatientData\\d*.csv", full.names=TRUE)
    if (is_empty(patientdatfiles)){
      stop("'patientdatfolder' contains no files with extension .csv OR 'patientdatfolder' does not point to a valid path")
    }

    if (length(patientdatfiles) != length(stage2files)){
      stop("patientdatfolder and stage2folder contain different numbers of csv files")
    }

    col_specs =
      cols(.default = col_double(),
           estMTDCode = col_character(),
           RP2DCode = col_character())

    trial_summary <-
      map(stage2files,
          read_csv,
          col_types = col_specs, progress = FALSE) %>%
      bind_rows()

    patient_summary <-
      map(patientdatfiles,
          read_csv,
          col_types = col_specs, progress = FALSE) %>%
      bind_rows()


    for (i in 1:length(dose_outcome_curves_list)) {
      newrow =
        tibble(Scenario = dose_outcome_curves_list[[i]]$scenario,
               `True DLT Probability` =  paste0("{",paste0(dose_outcome_curves_list[[i]]$tox_curve,collapse=","),"}"),
               `True Efficacy Probability` = paste0("{",paste0(dose_outcome_curves_list[[i]]$eff_curve,collapse=","),"}"),
               `True MTD` =
                 max(c(which(dose_outcome_curves_list[[i]]$tox_curve <= primary_objectives[["tox_target"]] +
                               primary_objectives[["tox_delta_no_exceed"]]),0)));

      curr_admiss = (dose_outcome_curves_list[[i]]$eff_curve >= primary_objectives[["eff_target"]]) &
        (dose_outcome_curves_list[[i]]$tox_curve <=  primary_objectives[["tox_target"]] +
           primary_objectives[["tox_delta_no_exceed"]]);
      if(sum(curr_admiss) == 0) {
        newrow <-
          mutate(newrow, `Acceptable/Desirable Dose` = as.character(0))
        curr_admiss = c(T, curr_admiss);
      } else if(sum(curr_admiss) == 1) {
        newrow <-
          mutate(newrow, `Acceptable/Desirable Dose` = as.character(which(curr_admiss)))
        curr_admiss = c(F, curr_admiss);
      } else {
        newrow <-
          mutate(newrow, `Acceptable/Desirable Dose` = paste0("{",paste0(which(curr_admiss),collapse=","),"}"))
        curr_admiss = c(F, curr_admiss);
      }

      generating_params_for_display <-
        bind_rows(generating_params_for_display, newrow);

      mtd_as_logical = (0:length(dose_outcome_curves_list[[i]]$tox_curve)) == max(c(which(dose_outcome_curves_list[[i]]$tox_curve <= primary_objectives[["tox_target"]] +
                                                                                       primary_objectives[["tox_delta_no_exceed"]]),0))
      generating_params =
        bind_rows(generating_params,
                  tibble(
                    scenario = dose_outcome_curves_list[[i]]$scenario,
                    dose_num = 0:length(dose_outcome_curves_list[[i]]$tox_curve),
                    true_dlt_prob = c(0, dose_outcome_curves_list[[i]]$tox_curve),
                    true_eff_prob = c(0, dose_outcome_curves_list[[i]]$eff_curve),
                    tox_target = primary_objectives["tox_target"] + primary_objectives["tox_delta_no_exceed"],
                    eff_target = primary_objectives["eff_target"],
                    is_mtd = mtd_as_logical,
                    is_acceptable = curr_admiss))
      scen <- c(scen, dose_outcome_curves_list[[i]]$scenario)
    }
  }

  ndose <- length(unique(generating_params[,"dose_num"]))-1

  #generating_params <- generating_params[order(generating_params[, "scenario"]),]
  #scen_num <- vector(length=nrow(generating_params))
  #for (j in 1:length(unique(generating_params[,"scenario"]))){
  #  myind <- which(generating_params[, "scenario"]==unique(generating_params[, "scenario"])[j])
  #  scen_num[myind] <- j
  #}
  #scen_designation <- ceiling(scen_num/scen_per_page)

  #generating_params <- cbind(generating_params, scen_designation)
  generating_params <-
    generating_params %>%
    arrange(scenario) %>%
    mutate(
      scennum =
        factor(scenario) %>%
        fct_inorder() %>%
        as.numeric(),
      scen_designation =
        ceiling(scennum / scen_per_page),
      scenario =
        factor(scenario) %>%
        fct_inorder() %>%
        fct_relabel(~paste0("Scenario ", .x))) %>%
    arrange(scenario, dose_num);

  generating_params_tall =
    generating_params %>%
    pivot_longer(true_dlt_prob:true_eff_prob, names_to = "type", values_to = "prob") %>%
    mutate(type =
             factor(type,
                    levels = c("true_dlt_prob","true_eff_prob"),
                    labels = c("DLT","Response")),
           is_acceptable_by_dose_num =
             factor((is_acceptable * (dose_num == 0)) +
                      2 * ((1 - is_acceptable) * (dose_num == 0)) +
                      3 * ((1 - is_acceptable) * (dose_num > 0)) +
                      4 * (is_acceptable * (dose_num > 0)),
                    levels = c(1,2,3,4))) %>%
    arrange(scenario, dose_num);

  if(is.null(design_labels)){
    design_labels <- c(1:length(unique(trial_summary$design)))
  } else {
    if (!near(length(design_labels), length(unique(trial_summary$design)))){
      stop("'design_labels' must be the same length as the number of unique designs")
    }
  }

  # Elizabeth: can you confirm that the following block of code is equivalent
  # to the commented code that immediately follows?
  trial_summary <-
    arrange(trial_summary, design, scenario) %>%
    mutate(designnum =
             factor(design) %>%
             fct_inorder() %>%
             as.numeric(),
           set_designation = ceiling(designnum / design_per_page),
           scennum =
             factor(scenario) %>%
             fct_inorder() %>%
             as.numeric(),
           scen_designation = ceiling(scennum / scen_per_page))
  #for (j in 1:length(unique(trial_summary$design))){
  #  trial_summary$designnum[trial_summary$design==unique(trial_summary$design)[j]] <- j
  #}
  #trial_summary$set_designation <- ceiling(trial_summary$designnum/design_per_page)

  #trial_summary <- arrange(trial_summary, scenario)
  #for (j in 1:length(unique(trial_summary$scenario))){
  #  trial_summary$scennum[trial_summary$scenario==unique(trial_summary$scenario)[j]] <- j
  #}
  #trial_summary$scen_designation <- ceiling(trial_summary$scennum/scen_per_page)

  trial_summary =
    trial_summary %>%
    arrange(design, scenario, array_id, sim_id) %>%
    mutate(design_label =
             factor(design,
                    levels = unique(design),
                    labels = design_labels,
                    ordered = T),
           scenario =
             factor(scenario) %>%
             fct_inorder() %>%
             fct_relabel(~paste0("Scenario ", .x)))

  # Result 1: Distribution of potential recommended dose levels (coarsened)

  trial_summary_RP2D =
    trial_summary %>%
    mutate(RP2DCode_truth =
             case_when(
               RP2DCode != "2Y"  & RP2DAcceptable == 1 ~ "NR",
               RP2DCode != "2Y"  & RP2DAcceptable == 0 ~ "NW",
               RP2DCode == "2Y" & RP2DAcceptable == 0 ~ "RW",
               TRUE ~ "RR"
             )) %>%
    group_by(set_designation, scen_designation, design, scenario,  design_label, bestP2D, RP2DCode_truth) %>%
    tally %>%
    mutate(sum_n = sum(n)) %>%
    ungroup() %>%
    mutate(prop_n = n / sum_n) %>%
    mutate(RP2DAcceptable = factor(RP2DCode_truth %in% c("NR","RR")),
           RP2DCode = factor(RP2DCode_truth,
                             levels = c("NR","NW","RW","RR"),
                             labels = c("No Rec\n(correct)",
                                        "No Rec\n(wrong)",
                                        "Rec\n(unaccept)",
                                        "Rec\n(accept)"),
                             ordered = T)) %>%
    arrange(set_designation, scen_designation, design, scenario,  RP2DCode) %>%
    group_by(set_designation, scen_designation, design, scenario,  design_label) %>%
    mutate(text_height =
             1 - (c(0,cumsum(prop_n)[-length(prop_n)]) + prop_n/2)) %>%
    ungroup()

  gen_param_plot <-
    vector("list", length = length(unique(generating_params_tall$scen_designation)))

  for (j in unique(generating_params_tall$scen_designation)){
    subdat <- filter(generating_params_tall, scen_designation==j)
    if (setequal(c("1", "3"), subdat$is_acceptable_by_dose_num)){
      outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(4,1)];
    } else if (setequal(c("2", "4"), subdat$is_acceptable_by_dose_num)){
      outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(2,5)];
    } else if (setequal(c("2", "3", "4"), subdat$is_acceptable_by_dose_num)){
      outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(2,1,5)];
    } else if (setequal(c("1", "2", "3", "4"), subdat$is_acceptable_by_dose_num)){
      outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(4,2,1,5)];
    }
    gen_param_plot[[j]] <-
      ggplot(subdat,
             aes(x = dose_num)) +
      geom_point(aes(y = prob, col = is_acceptable_by_dose_num, shape = type),
                 size = 3) +
      geom_line(data = filter(subdat, dose_num > 0),
                aes(y = prob, group = type), alpha = 0.25) +
      geom_point(data = filter(subdat, is_acceptable == 1),
                 aes(y = prob),
                 fill = "#00000000",
                 shape = 22,
                 size = 6) +
      facet_grid(scenario ~ ., scales = "free_y") +
      labs(x = "Dose Number",
           y = "",
           color = "Acceptable\nDose",
           linetype = "Endpoint",
           shape = "Endpoint") +
      scale_y_continuous(expand = expansion(mult = 0.05)) +
      scale_color_manual(values = outcome_colors) +
      guides(col = FALSE,
             shape = guide_legend(nrow = 1)) +
      theme(text = element_text(size = 10),
            legend.position = "top",
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_text_size),
            legend.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "pt"),
            legend.spacing = unit(1,units = "pt"),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank());
  }

  acc_dose_rec_plot <- vector("list", length = 0)
  for (j in unique(trial_summary_RP2D$set_designation)){
    for (k in unique(trial_summary_RP2D$scen_designation)){
      subdat <- filter(trial_summary_RP2D, set_designation==j & scen_designation==k)
      if (setequal(c("No Rec\n(correct)", "Rec\n(unaccept)"), subdat$RP2DCode)){
        outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(4,1)];
      } else if (setequal(c("No Rec\n(wrong)","Rec\n(accept)"), subdat$RP2DCode)){
        outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(2,5)];
      } else if (setequal(c("No Rec\n(wrong)","Rec\n(unaccept)","Rec\n(accept)"), subdat$RP2DCode)){
        outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(2,1,5)];
      } else if (setequal(c("No Rec\n(correct)","No Rec\n(wrong)","Rec\n(unaccept)","Rec\n(accept)"), subdat$RP2DCode)){
        outcome_colors = RColorBrewer::brewer.pal(5, "RdYlGn")[c(4,2,1,5)];
      }
      myplot =
        ggplot(data = subdat,
               aes(x = 1,
                   y = prop_n,
                   group = interaction(RP2DAcceptable, RP2DCode))) +
        geom_bar(aes(fill = RP2DCode),
                 stat = "identity",
                 color = NA) +
        geom_text(data = filter(subdat, prop_n > min_prop_to_write),
                  aes(x = 1,
                      y = text_height,
                      label = paste0(formatC(100 * prop_n, digits = 1, format = "f"),"%")),
                  size = text_size) +
        geom_col(aes(color = RP2DAcceptable,
                     size = RP2DAcceptable),
                 fill = "#FFFFFF00") +
        facet_grid(scenario ~ design_label, scales = "free_y",switch="both") +
        scale_y_reverse(labels = NULL, expand = expansion(add = 0.01)) +
        scale_x_continuous(expand = expansion(add = 0.002)) +
        scale_fill_manual(values = outcome_colors) +
        scale_color_manual(values = c("#FFFFFF00","black"), labels = c("No", "Yes")) +
        scale_size_manual(values = c(0.5, 0.75), labels = c("No", "Yes")) +
        labs(x="Design",
             y="",
             color = "Good\nOutcome",
             size = "Good\nOutcome",
             fill = "Outcome") +
        guides(color = guide_legend(nrow = 1),
               fill = guide_legend(nrow = 1)) +
        theme(text = element_text(size = 10),
              legend.position = "top",
              legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_text_size),
              legend.margin = margin(t = 0, r = 8, b = 0, l = 0, unit = "pt"),
              legend.spacing = unit(1,units = "pt"),
              strip.text.x = element_text(margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "pt")),
              axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title.y = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank());
      subplot <- list(myplot)
      acc_dose_rec_plot <- c(acc_dose_rec_plot, subplot)
      rm(list=c("myplot", "subplot"))
    }
  }

  dose_rec <-
    filter(trial_summary_RP2D, RP2DAcceptable == TRUE) %>%
    select(scenario, design_label, prop_n)
  dose_rec_table <-
    dose_rec %>%
    pivot_wider(id_cols = scenario,
                names_from = design_label,
                values_from = prop_n)

  samp_table <-
    trial_summary %>%
    group_by(scenario,  design_label) %>%
    summarize(mean_patients = mean(n_total_enrolled)) %>%
    pivot_wider(id_cols = scenario,
                names_from = design_label,
                values_from = mean_patients)

  patient_summary <-
    arrange(patient_summary, design, scenario) %>%
    mutate(dose_rec =
             (dose_num == RP2D),
           design_label =
             factor(design,
                    levels = unique(design),
                    labels = design_labels,
                    ordered = T),
           scenario =
             factor(scenario,
                    levels = unique(scenario)[order(unique(scenario))],
                    labels = paste0("Scenario ", unique(scenario)[order(unique(scenario))])),
           designnum =
             factor(design) %>%
             fct_inorder() %>%
             as.numeric(),
           set_designation = ceiling(designnum / design_per_page),
           scennum =
             factor(scenario) %>%
             fct_inorder() %>%
             as.numeric(),
           scen_designation = ceiling(scennum / scen_per_page))


  num_at_RP2D <-
    patient_summary %>%
    group_by(scenario, design_label, array_id, sim_id) %>%
    summarize(atRP2D = sum(dose_rec))

  # Result 2: Distribution of dose levels at 10th, 30th, 50th patient

  patient_summary_num_atDoseLevels =
    patient_summary %>%
    group_by(design, scenario,  design_label, subj_id, set_designation, scen_designation, dose_num) %>%
    tally %>%
    mutate(sum_n = sum(n)) %>%
    ungroup() %>%
    mutate(prop_n = n / sum_n) %>%
    arrange(design, scenario, subj_id, dose_num) %>%
    group_by(design, scenario,  design_label, subj_id, set_designation, scen_designation) %>%
    mutate(text_height = 1 - (c(0,cumsum(prop_n)[-length(prop_n)]) + prop_n/2)) %>%
    ungroup() %>%
    left_join(y = generating_params) %>%
    mutate(dose_num = factor(dose_num,
                             ordered = T),
           is_acceptable = factor(is_acceptable,
                                  levels = c(0, 1),
                                  labels = c("No","Yes")));


  dose_over_time_plot <- vector("list", length = 0)
  for (j in unique(patient_summary_num_atDoseLevels$set_designation)) {
    for (k in unique(patient_summary_num_atDoseLevels$scen_designation)) {

      subdat <-
        filter(patient_summary_num_atDoseLevels,
               set_designation == j,
               scen_designation == k,
               subj_id %in% c(10, 30, 50))

      myplot =
        ggplot(data = subdat,
               aes(x = factor(subj_id))) +
        geom_col(aes(y = prop_n,
                     group = interaction(is_acceptable, dose_num),
                     fill = dose_num),
                 color = NA) +
        geom_col(aes(y = prop_n,
                     group = interaction(is_acceptable, dose_num),
                     color = is_acceptable,
                     size = is_acceptable),
                 fill = "#FFFFFF00") +
        geom_text(data = filter(subdat, prop_n > min_prop_to_write),
                  aes(y = text_height,
                      label = paste0(formatC(100 * prop_n, digits = 1, format = "f"),"%")),
                  size = text_size) +
        facet_grid(scenario ~ design_label, scales = "free_y", switch = "y") +
        scale_y_reverse(labels = NULL, expand = expansion(add = 0.01)) +
        scale_x_discrete(expand = expansion(add = 0.002)) +
        scale_fill_brewer(palette = "Blues") +
        scale_color_manual(values = c("#FFFFFF00","black")) +
        scale_size_manual(values = c(0.5, 0.75)) +
        labs(x="Subj ID",
             y="",
             color = "Acceptable",
             size = "Acceptable",
             fill = "Dose Number") +
        guides(color = guide_legend(nrow = 1),
               size = guide_legend(nrow = 1),
               fill = guide_legend(nrow = 1)) +
        theme(text = element_text(size = 10),
              legend.position = "top",
              legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_text_size),
              legend.margin = margin(t = 0, r = 3, b = 0, l = 0, unit = "pt"),
              legend.spacing = unit(1,units = "pt"),
              strip.text.x = element_text(margin = margin(t = 3, r = 0, b = 3, l = 0, unit = "pt")),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank());
      subplot <- list(myplot)
      dose_over_time_plot <- c(dose_over_time_plot, subplot)
      rm(list=c("myplot", "subplot"))
    }
  }

  prop_accept <-
    left_join(patient_summary, generating_params, by=c("scenario", "dose_num")) %>%
    group_by(scenario, design_label, array_id, sim_id) %>%
    summarize(total = n(), ataccept = sum(is_acceptable)) %>%
    mutate(prop_acc = ataccept/total)

  samp_plot <- ggplot(data=trial_summary, aes(x=scenario, y = n_total_enrolled, fill=design_label)) +
    geom_boxplot(color="black", lwd=0.2, position = position_dodge2(padding=0.3)) +
    ylab("Number of Patients Enrolled") + xlab("Scenario") + scale_fill_discrete(name="Design")

  samp_RP2D_plot <- ggplot(data=num_at_RP2D, aes(x=scenario, y = atRP2D, fill=design_label)) +
    geom_boxplot(color="black", lwd=0.2, position = position_dodge2(padding=0.3)) +
    ylab("Number of Patients at RP2D") + xlab("Scenario") + scale_fill_discrete(name="Design")

  prop_accep_plot <- ggplot(data=prop_accept, aes(x=scenario, y = prop_acc, fill=design_label)) +
    geom_boxplot(color="black", lwd=0.2, position = position_dodge2(padding=0.3)) +
    ylab("Proportion of Patients \n Receiving Acceptable Dose") + xlab("Scenario") + scale_fill_discrete(name="Design")

  tables <- list(
    gen_params_table = generating_params_for_display,
    acc_dose_rec_table = dose_rec_table,
    mean_patients_table = samp_table
  )

  plots <- list(
    gen_params_plot = gen_param_plot, #gen_param_plot[[1]],
    acc_dose_rec_plot = acc_dose_rec_plot, #acc_dose_rec_plot[[1]],
    dose_over_time_plot = dose_over_time_plot,
    n_patients_plot = samp_plot,
    n_patients_RP2D_plot = samp_RP2D_plot,
    prop_patients_accep_plot = prop_accep_plot
  )

  results <- list(
    tables = tables,
    plots = plots
  )

  return(results)
}
