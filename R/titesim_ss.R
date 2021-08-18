#' Expanded version of the dfcrm::titesim function to incorporate design elements into the
#' time-to-event continual reassessment method that are especially useful for
#' seamless simulations.
#'
#' This is the simulator function for a TITE-CRM trial. It is meant to be called in the context of
#' twostage_simulator() rather than by the user directly. Required input includes both the design elements as well as the
#' true, dose-toxicity curve that is
#' generally unknown in the real world. The simulator runs a certain number of
#' simulated trials, creating data according to the dose-toxicity curve and
#' making assignments according to the titecrm model. Various operating characteristics
#' are reported.
#'
#' @param PI Numeric vector with entries between 0 and 1 of true toxicity probabilities; assumed to match
#' the order of the dose labels.
#' @param prior Numeric vector with entries between 0 and 1 of anticipated toxicity probabilities, assumed to match
#' the order of the dose labels. More commonly and accurately called the skeleton.
#' @param target Scalar value between 0 and 1 giving the targeted rate of DLT.
#' @param n A positive integer indicating the maximum number of patients to enroll.
#' @param x0 A positive integer indicating the starting dose level.
#' @param nsim A positive integer indicating the number of simulated trials to conduct.
#' @param restrict A logical value. If TRUE, indicates that safety constraints should be enacted. There are four:
#' (1) No skipping doses in escalation; (2) No escalation before followup of followup_b4_esc on at least one
#' patient at current or larger dose; (3) No assignment to dose with estimated DLT rate beyond no.exceed +
#' target; (4) Stopping trial altogether if at any point after patient earliest_stop the estimated DLT rate
#' at all dose levels exceeds no.exceed + target.
#' @param obswin A positive number indicating the number of units of time over which DLTs are defined.
#' @param tgrp This argument is an artefact of the dfcrm package and should be left blank.
#' @param rate A positive number indicating the number of patients expected per unit obswin.
#' @param accrual A string, either "fixed" or "poisson", denoting the patient accrual process.
#' @param surv A string, either "uniform" or "exponential", denoting the time-to-DLT distribution.
#' @param surv_rate A positive number. If surv is "exponential", this is the rate for the time-to-DLT
#' distribution.
#' @param scheme A string indicating the weighting scheme for patients who are free of DLT but have not
#' completed followup. Can be either "polynomial", "logistic", or "adaptive."
#' @param scheme_args A named list with elements "scheme power" (if scheme = "polynomial"), "scheme_int"
#' and "scheme_slope" (if scheme = "logistic"), or no elements (if scheme = "adaptive").
#' @param count A logical value; if TRUE, the progress of the simulations will be plotted.
#' @param method A string indicating the method for fitting model. The original titesim function allows
#' "mle" or "bayes"; here, only "bayes" is allowed.
#' @param model A string indicating the type of model. The original titesim function allows "empiric"
#' (sometimes known as the power model) or "logistic".
#' @param intcpt A numeric value giving the intercept parameter when using the "logistic" model.
#' @param scale A numeric value giving the prior standard deviation on the parameter beta.
#' @param seed A positive integer random seed.
#' @param conf.level A number between 0 and 1 indicating the confidence limits to report.
#' @param no.exceed A positive number indicating by how much the toxicity rate can exceed the DLT threshold
#' before the dose is unacceptable. To be more clear, no dose will be assigned with estimated pr(DLT) >
#' target + no.exceed, and the trial will stop altogether if this holds for dose level 1.
#' @param cohort.size A positive integer indicating the size of cohorts between successive queries of the
#' CRM model.
#' @param first.cohort.only A logical value indicating if cohort.size applies only to the first cohort.
#' @param n.at.MTD A positive integer which is the number of patients that should be at the estimated MTD
#' before the trial stops. If n.at.MTD patients have been assigned to a single dose, and this dose is the
#' current MTD that would be assigned to the next patient, then the trial stops.
#' @param followup_b4_esc A positive number indicating how much follow-up time at least one patient has to have
#' on a given dose before the model escalates to the next dose.
#' @param earliest_stop A positive integer indicating how many patients must be observed before the trial can stop.
#' @return The function returns a list, with named components last_sim and all_sim. Access the last_sim to
#' see the individual details of the single last simulation, or access all_sim to see summaries across all
#' simulations performed.
#' @references
#'
#' \insertRef{boonstra2020}{seamlesssim}
#'
#' \insertRef{dfcrm2019}{seamlesssim}
#'
#' @importFrom stats rexp binomial runif
#' @import dfcrm
#' @export
titesim_ss = function (PI, prior, target, n, x0, nsim = 1, restrict = TRUE,
                         obswin = 1, tgrp = obswin, rate = 1, accrual = "fixed", surv = "uniform", surv_rate = obswin/10,
                         scheme = "polynomial", scheme_args = list(scheme_power = 1), count = TRUE, method = "bayes", model = "empiric",
                         intcpt = 3, scale = sqrt(1.34), seed = 1009,
                         conf.level = 0.5, no.exceed = Inf, cohort.size = 1, first.cohort.only = T, n.at.MTD = Inf, followup_b4_esc = obswin, earliest_stop = 6)
{
  set.seed(seed);
  #number exposed, number toxicities
  nexpt <- ntox <- matrix(0,nsim,length(prior));
  #which dose selected
  sel <- matrix(0, nsim, length(prior)+1);
  #parameter estimates
  BETAHAT <- matrix(rep(NA, nsim * n), nrow = nsim)
  #final parameter estimate, duration of trial, number enrolled
  final.est <- DURATION <- n.enrolled <- rep(NA, nsim)
  ###Begin Phil's modification
  #after which patient was trial stopped for toxicity? (0 if never stopped)
  stop.for.tox = numeric(nsim);
  #size of initial cohort
  first.cohort.size = cohort.size;
  ###End Phil's modification
  for (r in 1:nsim) {
    cohort.size = first.cohort.size;###Phil's modification
    if (count) {
      cat("simulation number:", r, "\n")
    }
    if (accrual == "fixed") {
      next_arrival <- obswin/rate
    } else if (accrual == "poisson") {
      next_arrival <- rexp(1, rate/obswin)
    }
    if (length(x0) > 1) {
      stop("Phil's modifications only made for standard 1-stage TITE-CRM")
    } else {
      if (method == "mle") {
        stop(" Require an initial design for mle-CRM!")
      }
      bethat <- 0;
      time_tox <- all_tox <- level <- arrival <- numeric(n);
      cur <- x0;
      i=1;
      while(i<=n) {
        arrival[i] <- next_arrival;
        level[i] <- cur;
        next_tox <- rbinom(1, 1, PI[cur])
        if (surv == "uniform") {
          time_next_tox <- runif(1, 0, obswin)/next_tox;
        } else if (surv == "exponential") {
          time_next_tox <- min(obswin,rexp(1,1/surv_rate))/next_tox;
        }
        all_tox[i] <- next_tox;
        time_tox[i] <- time_next_tox;
        trial_time_tox <- (time_tox + arrival)[1:i];

        if(i<n) {
          if (accrual == "fixed") {
            next_arrival <- next_arrival + obswin/rate
          } else if (accrual == "poisson") {
            next_arrival <- next_arrival + rexp(1, rate/obswin)
          }
          tox_obs <- rep(0, i)
          tox_obs[trial_time_tox <= next_arrival] <- 1
          followup <- pmin(pmin(next_arrival,trial_time_tox) - arrival[1:i], obswin)
          obj <- titecrm_ss(prior, target, tox_obs, level[1:i], followup = followup,
                              obswin = obswin, scheme = scheme, scheme_args = scheme_args,
                              method = method,
                              model = model, intcpt = intcpt, scale = scale,
                              var.est = T, conf.level = conf.level);
        } else {
          followup <- pmin(trial_time_tox-arrival,obswin)
          obj <- titecrm_ss(prior, target, all_tox, level, weights = rep(1,n), method = method, model = model, intcpt = intcpt,
                              scale = scale, var.est = T, conf.level = conf.level);
        }

        ###Begin seamless sim modifications
        if(i%%cohort.size==0) {#Only modify the dose at the end of each cohort
          if(restrict) {
            max.possible = max(c(0,which((obj$ptox-target)<=no.exceed)));#To indicate whether to stop early due to toxicity
            if(max.possible==0){
              if(i>=earliest_stop) {
                stop.for.tox[r] = i;#Variable to indicate at what patient the trial stopped
                bethat <- c(bethat, rep(-Inf,n-i));
                cur = 0;
                break;
              } else {
                max.possible = 1;#Continue the trial at the lowest dose if fewer than 6 patients have been enrolled
              }
            }
            if((!any(followup[which(level>=level[i])]>=followup_b4_esc))) {
              #Do not escalate until you've followed at least one patient for a certain period at this or a higher dose level
              max.possible = min(max.possible,cur)
            }
            cur <- min(obj$mtd, (cur + 1), max.possible)#Do not skip a dose, do not exceed all of the previous restrictions
          } else {
            cur <- obj$mtd
          }
        }

        if(i == n) {
          break;
        }

        if(i%%cohort.size == 0 & first.cohort.only == T) {cohort.size=1;}#If cohorts only apply to the first batch of patients, set subsequent cohort sizes to 1.


        if(n.at.MTD < Inf && sum(level==cur) >= n.at.MTD && sum((followup*(1-tox_obs)/obswin + tox_obs)[level[1:i]==cur][1:n.at.MTD]) >= n.at.MTD)  {
          obj <- titecrm_ss(prior, target, all_tox[1:i], level[1:i], weights = rep(1,i), method = method, model = model, intcpt = intcpt,
                              scale = scale, var.est = T, conf.level = conf.level);

          bethat = c(bethat, obj$est, rep(Inf,n-i-1));
          break;
        }
        ###End seamless sim modifications
        i=i+1;
        bethat <- c(bethat, obj$est);
      }
      BETAHAT[r, ] <- bethat;
      est <- obj$est;
      msg <- "Okay"
    }

    time_tox <- time_tox[1:i];
    all_tox <- all_tox[1:i];
    level <- level[1:i] ;
    arrival <- arrival[1:i];

    sel[r,cur+1] <- 1;#first index of 'sel' corresponds to stopping for toxicity
    final.est[r] <- est;
    DURATION[r] <- max(arrival) + obswin;
    n.enrolled[r] = i;
    for (k in 1:length(prior)) {
      nexpt[r,k] <- length(which(level == k))
      ntox[r,k] <- length(which(all_tox == 1 & level == k))
    }
  }
  if (length(x0) == 1) {
    design <- paste("TITE-CRM starting at dose", x0);
  } else {
    design <- "Two-stage TITE-CRM";
  }
  foo <- list(all_sim = list(PI = PI, prior = prior, target = target, n = n, n.enrolled = n.enrolled,
                             x0 = x0, nsim = nsim, MTD = colMeans(sel), level = colMeans(nexpt), tox = colMeans(ntox),
                             beta.hat = BETAHAT, final.est = final.est, Duration = DURATION,
                             design = design, method = method, prior.var = scale^2,
                             model = model, intcpt = intcpt, restriction = restrict,
                             seed = seed, tite = TRUE, dosescaled = obj$dosescaled,
                             msg = msg, obswin = obswin, tgrp = tgrp, rate = rate,
                             accrual = accrual, scheme = scheme, scheme_args = scheme_args,
                             no.exceed = no.exceed, cohort.size = first.cohort.size, first.cohort.only = first.cohort.only,
                             stop.for.tox = stop.for.tox, n.at.MTD = n.at.MTD, followup_b4_esc = followup_b4_esc, earliest_stop = earliest_stop,
                             sel = sel, nexpt = nexpt, ntox = ntox),
              last_sim = list(PI = PI, prior = prior, target = target, n = n, n.enrolled = n.enrolled,
                              x0 = x0, nsim = nsim, MTD = cur, level = level, tox = all_tox,
                              beta.hat = bethat, final.est = est, arrival = arrival,
                              toxicity.time = time_tox, toxicity.study.time = trial_time_tox, design = design,
                              method = method, prior.var = scale^2, model = model,
                              intcpt = intcpt, restriction = restrict, seed = seed,
                              tite = TRUE, dosescaled = obj$dosescaled, msg = msg,
                              obswin = obswin, tgrp = tgrp, rate = rate, accrual = accrual,
                              scheme = scheme, scheme_args = scheme_args,
                              post.var = obj$post.var, ptox = obj$ptox,
                              ptoxL = obj$ptoxL, ptoxU = obj$ptoxU, conf.level = obj$conf.level,
                              no.exceed = no.exceed, cohort.size = first.cohort.size, first.cohort.only = first.cohort.only,
                              stop.for.tox = stop.for.tox, n.at.MTD = n.at.MTD, followup_b4_esc = followup_b4_esc, earliest_stop = earliest_stop))

  class(foo) <- "titesim_ss"
  foo
}
