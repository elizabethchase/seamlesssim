#' Function to fit a TITE-CRM
#'
#' This is a function that makes dose recommendations for the next patient given the inputted data and
#' the design parameters according to the continual reasssessment method. At the end of the trial, it returns
#' estimates of the dose-toxicity curve under a one-parameter model where dose is the only predictor
#' and prints out a dose recommendation. IMPORTANT: NO SAFETY CONSTRAINTS ARE IMPLEMENTED
#' IN THIS FUNCTION. IT ONLY PRINTS OUT THE MODEL-BASED DOSE ASSIGNMENT FOR THE NEXT PATIENT;
#' IT IS UP TO THE USER TO  DETERMINE WHETHER ALL SAFETY CONSTRAINTS WOULD BE SATISFIED BY ANY
#' GIVEN DOSE ASSIGNMENT AND TO REDUCE THE ASSIGNMENT AS NECESSARY. This function is based upon
#' the titecrm function by Ken Cheung.
#'
#' @param prior A numeric vector with values between 0 and 1; the anticipated probabilities
#' of toxicity for each dose. More commonly called the skeleton.
#' @param target A scalar between 0 and 1 giving the targeted rate of DLT.
#' @param tox An integer vector of 0s and 1s the same length as the current number of patients enrolled,
#' indicating whether or not that patient had a toxicity.
#' @param level An integer vector of dose numbers indicating dose assignments for all currently
#' enrolled patients. Same length as tox.
#' @param n An integer greater than 0 indicating the number of patients already enrolled, equal
#' to the lengths of tox and level.
#' @param weights A numeric vector of weights between 0 and 1 that control the likelihood
#' contribution for each patient, in the situation where different patients are observed for different
#' lengths of time. Same length as tox.
#' @param followup A positive numeric vector indicating the number of units of time that each
#' patient has been followed; same length as tox.
#' @param entry Positive numeric vectors of entry and exit times; alternative to calculating followup.
#' Same length as tox.
#' @param exit Positive numeric vectors of entry and exit times; alternative to calculating followup.
#' Same length as tox.
#' @param obswin A positive numeric value indicating the number of units of time over which DLTs are defined.
#' @param scheme A string indicating the weighting scheme for patients who are free of DLT but
#' have not completed followup. Must be either "polynomial", "logistic", or "adaptive". "polynomial" is the
#' default.
#' @param scheme_args A named list with elements "scheme_power" (if "scheme" = "polynomial"),
#' "scheme_int" and "scheme_slope" (if "scheme" = "logistic"), or no elements (if "scheme" = "adaptive").
#' @param conf.level A number between 0 and 1; the confidence limits to report. Default is 0.9.
#' @param dosename A vector the same length as prior giving a list of names/identifiers for the different
#' doses.
#' @param include From titecrm documentation: "A subset of patients included in the dose calculation".
#' Default is to include all patients.
#' @param pid A vector of length n giving each patient's identifier. Default is to assign each patient
#' an identifier from 1 to n.
#' @param method A string indicating the method for fitting the model. The original titecrm
#' function allows "mle" or "bayes"; titecrm_phil only includes "bayes".
#' @param model A string indicating the type of model. The original titesim function allows
#' "empiric" (sometimes known as the power model) or "logistic"; titecrm_phil only includes "empiric".
#' @param var.est A logical value indicating if the posterior variance of model parameters should be returned.
#' Default is TRUE.
#' @param scale A positive numeric value indicating the prior standard deviation on the parameter beta.
#' Default is the square root of 1.34.
#' @param intcpt A fixed numeric value of the intercept parameter when using the "logistic" model.
#' Default is 3.
#' @param model.detail From titecrm documentation: "If FALSE, the model content of an mtd
#' object will not be displayed. Default is TRUE".
#' @param patient.detail From titecrm documentation: "If FALSE, patient summary of an mtd
#' object will not be displayed. Default is TRUE".
#' @param tite From titecrm documentation: "If FALSE, the time components in patient summary of
#' an mtd object will be omitted. Default is TRUE".
#' @return A named list with entries prior, target, tox, level, dosename, weights, followup, entry,
#' exit, obswin, scheme, scheme_args, model, method, model.detail, intcpt, conf.level, include,
#' tite, dosescaled, and patient.detail as described above, along with:
#' \describe{
#'   \item{prior.var}{The prior variance of beta, or the user-inputted value scale squared.}
#'   \item{post.var}{The posterior variance of beta.}
#'   \item{subset}{A vector of patient IDs indicating which patients were included in the dose calculation.}
#'   \item{estimate}{Posterior estimate of beta.}
#'   \item{mtd}{Estimated MTD.}
#'   \item{ptox}{Probability of toxicity at each dose.}
#'   \item{ptoxL}{Lower confidence interval bound on the probability of toxicity at each dose.}
#'   \item{ptoxU}{Upper confidence interval bound on the probability of toxicity at each dose.}
#'   \item{dosescaled}{Scaled dose levels.}
#' }
#' @references Ken Cheung, "dfcrm: Dose-Finding by the Continual Reassessment Method." Version 0.2-2.1, 2019.
#' https://CRAN.R-project.org/package=dfcrm
#' @importFrom stats optimize integrate qnorm
#' @import dfcrm
#' @export
titecrm_phil = function (prior, target, tox, level, n = length(level), weights = NULL,
                         followup = NULL, entry = NULL, exit = NULL, obswin = NULL,
                         scheme = "polynomial", scheme_args = list(scheme_power=1), conf.level = 0.9, dosename = NULL, include = 1:n,
                         pid = 1:n, method = "bayes", model = "empiric", var.est = TRUE,
                         scale = sqrt(1.34), intcpt = 3, model.detail = TRUE, patient.detail = TRUE,
                         tite = TRUE)
{
  if (is.null(weights)) {
    if (is.null(followup)) {
      followup <- exit - entry
    }
    if (scheme == "polynomial") {
      weights <- (followup^scheme_args$scheme_power)/(obswin^scheme_args$scheme_power);
    } else if (scheme == "logistic") {
      weights <- ifelse(followup>=obswin,1,1/(1+exp(-scheme_args$scheme_int - followup * scheme_args$scheme_slope)));
    } else if (scheme == "adaptive") {
      support <- sort(followup[tox == 1])
      z <- length(support)
      if (z) {
        for (i in 1:n) {
          m <- length(support[support <= followup[i]])
          if (!m)
            weights[i] <- followup[i]/support[1]/(z +
                                                    1)
          else if (m == z)
            weights[i] <- (z + (followup[i] - support[z])/(obswin -
                                                             support[z]))/(z + 1)
          else weights[i] <- (m + (followup[i] - support[m])/(support[m +
                                                                        1] - support[m]))/(z + 1)
        }
      } else {
        weights <- followup/obswin
      }
    } else {
      stop(" Weighting scheme undefined!")
    }
    weights <- pmin(weights, 1)
  }
  if (any(weights > 1) | any(weights < 0))
    stop(" Weights have to be between 0 and 1!")
  if (is.null(pid)) {
    if (!(length(tox) == length(level) & length(tox) == length(weights)))
      stop(" tox, level, and weights are of different lengths!")
  }
  else {
    if (!(length(tox) == length(level) & length(tox) == length(weights) &
          length(tox) == length(pid)))
      stop(" pid, tox, level, and weights are of different lengths!")
  }
  weights[tox == 1] <- 1;
  y1p <- tox[include]
  w1p <- weights[include]
  if (model == "empiric") {
    dosescaled <- prior
    x1p <- prior[level[include]]
    if (method == "mle") {
      if (sum(y1p) == 0 | sum(y1p) == length(y1p))
        stop(" mle does not exist!")
      est <- optimize(lcrm, c(-10, 10), x1p, y1p, w1p,
                      tol = 1e-04, maximum = TRUE)$max
      if (var.est) {
        e2 <- integrate(crmht2, -10, 10, x1p, y1p, w1p,
                        500, abs.tol = 0)[[1]]/integrate(crmh, -10,
                                                         10, x1p, y1p, w1p, 500, abs.tol = 0)[[1]]
      }
    }
    else if (method == "bayes") {
      den <- integrate(crmh, -Inf, Inf, x1p, y1p, w1p,
                       scale, abs.tol = 0)[[1]]
      est <- integrate(crmht, -10, 10, x1p, y1p, w1p, scale,
                       abs.tol = 0)[[1]]/den
      if (var.est) {
        e2 <- integrate(crmht2, -10, 10, x1p, y1p, w1p,
                        scale, abs.tol = 0)[[1]]/den
      }
    }
    else {
      stop(" unknown estimation method")
    }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2 - est^2
      crit <- qnorm(0.5 + conf.level/2)
      lb <- est - crit * sqrt(post.var)
      ub <- est + crit * sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model == "logistic") {
    dosescaled <- log(prior/(1 - prior)) - intcpt
    if (!all(dosescaled < 0)) {
      stop("intercept parameter in logit model is too small: scaled doses > 0!")
    }
    x1p <- dosescaled[level[include]]
    if (method == "mle") {
      if (sum(y1p) == 0 | sum(y1p) == length(y1p))
        stop(" mle does not exist!")
      est <- optimize(lcrmlgt, c(-10, 10), x1p, y1p, w1p,
                      intcpt, tol = 1e-04, maximum = TRUE)$max
      if (var.est) {
        e2 <- integrate(crmht2lgt, -10, 10, x1p, y1p,
                        w1p, 500, intcpt, abs.tol = 0)[[1]]/integrate(crmhlgt,
                                                                      -10, 10, x1p, y1p, w1p, 500, abs.tol = 0)[[1]]
      }
    }
    else if (method == "bayes") {
      den <- integrate(crmhlgt, -Inf, Inf, x1p, y1p, w1p,
                       scale, intcpt, abs.tol = 0)[[1]]
      est <- integrate(crmhtlgt, -10, 10, x1p, y1p, w1p,
                       scale, intcpt, abs.tol = 0)[[1]]/den
      if (var.est) {
        e2 <- integrate(crmht2lgt, -10, 10, x1p, y1p,
                        w1p, scale, intcpt, abs.tol = 0)[[1]]/den
      }
    }
    else {
      stop(" unknown estimation method")
    }
    ptox <- (1 + exp(-intcpt - exp(est) * dosescaled))^{
      -1
    }
    if (var.est) {
      post.var <- e2 - est^2
      crit <- qnorm(0.5 + conf.level/2)
      lb <- est - crit * sqrt(post.var)
      ub <- est + crit * sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt - exp(ub) * dosescaled))^{
        -1
      }
      ptoxU <- (1 + exp(-intcpt - exp(lb) * dosescaled))^{
        -1
      }
    }
  }
  else {
    stop(" model specified not available.")
  }
  if (all(ptox <= target)) {
    rec <- length(prior)
  }
  else if (all(ptox >= target)) {
    rec <- 1
  }
  else {
    rec <- order(abs(ptox - target))[1]
  }
  if (!var.est) {
    post.var <- ptoxL <- ptoxU <- NA
  }
  foo <- list(prior = prior, target = target, tox = tox, level = level,
              dosename = dosename, subset = pid[include], estimate = est,
              weights = weights, followup = followup, entry = entry,
              exit = exit, obswin = obswin, scheme = scheme, scheme_args = scheme_args, model = model,
              prior.var = scale^2, post.var = post.var, method = method,
              mtd = rec, include = include, pid = pid, model.detail = model.detail,
              intcpt = intcpt, ptox = ptox, ptoxL = ptoxL, ptoxU = ptoxU,
              conf.level = conf.level, patient.detail = patient.detail,
              tite = tite, dosescaled = dosescaled)
  class(foo) <- "mtd"
  foo
}

