% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/titecrm_ss.R
\name{titecrm_ss}
\alias{titecrm_ss}
\title{Expanded version of the dfcrm::titecrm function to incorporate some useful design elements of the
time-to-event continual reassessment method}
\usage{
titecrm_ss(
  prior,
  target,
  tox,
  level,
  n = length(level),
  weights = NULL,
  followup = NULL,
  entry = NULL,
  exit = NULL,
  obswin = NULL,
  scheme = "polynomial",
  scheme_args = list(scheme_power = 1),
  conf.level = 0.9,
  dosename = NULL,
  include = 1:n,
  pid = 1:n,
  method = "bayes",
  model = "empiric",
  var.est = TRUE,
  scale = sqrt(1.34),
  intcpt = 3,
  model.detail = TRUE,
  patient.detail = TRUE,
  tite = TRUE
)
}
\arguments{
\item{prior}{A numeric vector with values between 0 and 1; the anticipated probabilities
of toxicity for each dose. More commonly called the skeleton.}

\item{target}{A scalar between 0 and 1 giving the targeted rate of DLT.}

\item{tox}{An integer vector of 0s and 1s the same length as the current number of patients enrolled,
indicating whether or not that patient had a toxicity.}

\item{level}{An integer vector of dose numbers indicating dose assignments for all currently
enrolled patients. Same length as tox.}

\item{n}{An integer greater than 0 indicating the number of patients already enrolled, equal
to the lengths of tox and level.}

\item{weights}{A numeric vector of weights between 0 and 1 that control the likelihood
contribution for each patient, in the situation where different patients are observed for different
lengths of time. Same length as tox.}

\item{followup}{A positive numeric vector indicating the number of units of time that each
patient has been followed; same length as tox.}

\item{entry}{Positive numeric vectors of entry and exit times; alternative to calculating followup.
Same length as tox.}

\item{exit}{Positive numeric vectors of entry and exit times; alternative to calculating followup.
Same length as tox.}

\item{obswin}{A positive numeric value indicating the number of units of time over which DLTs are defined.}

\item{scheme}{A string indicating the weighting scheme for patients who are free of DLT but
have not completed followup. Must be either "polynomial", "logistic", or "adaptive". "polynomial" is the
default.}

\item{scheme_args}{A named list with elements "scheme_power" (if "scheme" = "polynomial"),
"scheme_int" and "scheme_slope" (if "scheme" = "logistic"), or no elements (if "scheme" = "adaptive").}

\item{conf.level}{A number between 0 and 1; the confidence limits to report. Default is 0.9.}

\item{dosename}{A vector the same length as prior giving a list of names/identifiers for the different
doses.}

\item{include}{From titecrm documentation: "A subset of patients included in the dose calculation".
Default is to include all patients.}

\item{pid}{A vector of length n giving each patient's identifier. Default is to assign each patient
an identifier from 1 to n.}

\item{method}{A string indicating the method for fitting the model. The original titecrm
function allows "mle" or "bayes". The choice of "mle" has not been tested extensively in titecrm_ss}

\item{model}{A string indicating the type of model. The original titesim function allows
"empiric" (sometimes known as the power model) or "logistic".  The choice of "logistic"
has not been tested extensively in titecrm_ss}

\item{var.est}{A logical value indicating if the posterior variance of model parameters should be returned.
Default is TRUE.}

\item{scale}{A positive numeric value indicating the prior standard deviation on the parameter beta.
Default is the square root of 1.34.}

\item{intcpt}{A fixed numeric value of the intercept parameter when using the "logistic" model.
Default is 3.}

\item{model.detail}{From titecrm documentation: "If FALSE, the model content of an mtd
object will not be displayed. Default is TRUE".}

\item{patient.detail}{From titecrm documentation: "If FALSE, patient summary of an mtd
object will not be displayed. Default is TRUE".}

\item{tite}{From titecrm documentation: "If FALSE, the time components in patient summary of
an mtd object will be omitted. Default is TRUE".}
}
\value{
A named list with entries prior, target, tox, level, dosename, weights, followup, entry,
exit, obswin, scheme, scheme_args, model, method, model.detail, intcpt, conf.level, include,
tite, dosescaled, and patient.detail as described above, along with:
\describe{
  \item{prior.var}{The prior variance of beta, or the user-inputted value scale squared.}
  \item{post.var}{The posterior variance of beta.}
  \item{subset}{A vector of patient IDs indicating which patients were included in the dose calculation.}
  \item{estimate}{Posterior estimate of beta.}
  \item{mtd}{Estimated MTD.}
  \item{ptox}{Probability of toxicity at each dose.}
  \item{ptoxL}{Lower confidence interval bound on the probability of toxicity at each dose.}
  \item{ptoxU}{Upper confidence interval bound on the probability of toxicity at each dose.}
  \item{dosescaled}{Scaled dose levels.}
}
}
\description{
This is a function that makes dose recommendations for the next patient given the inputted data and
the design parameters according to the continual reasssessment method. It is meant to be called in the context of
twostage_simulator() (by way of titesim_ss()) rather than by the user directly. At the end of the trial, it returns
estimates of the dose-toxicity curve under a one-parameter model where dose is the only predictor
and prints out a dose recommendation. IMPORTANT: No safety constraints are
implemented in this function. It only prints out the model-based dose assignment
for the next patient. The safety constraints are taken care of in
titesim_ss(). If the user is calling titecrm_ss directly, it is up to them to
determine whether all intended safety constraints are satisfied and to reduce the assignment as necessary.
}
\references{
\insertRef{boonstra2020}{seamlesssim}

\insertRef{dfcrm2019}{seamlesssim}
}
