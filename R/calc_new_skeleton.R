#' Function to update a CRM skeleton in a seamless design where the second
#' stage uses the CRM.
#'
#' This function updates an initial CRM skeleton when using the empiric (power) model, after
#' some data have been gathered on the dose-toxicity curve. This code is currently not used,
#' but is included for posterity.
#'
#' @param old_skeleton A numeric vector with entries between 0 and 1; the initial skeleton.
#' @param tox A vector with binary entries {0,1} representing the outcomes of the patients;
#' should be the same length as the number of observations collected.
#' @param level A positive integer vector of labels with length equal to that of tox indicating
#' the dose level that each observation was assigned to. Every element of this vector should,
#' in addition to being positive, be no larger than the length of skeleton, i.e. the number of
#' dose levels.
#' @param offset_seq A positive numeric vector of arbitrary length that should probably have
#' some values above and below 1. The old_skeleton corresponds to a value of 1, and including
#' offset_seq helps to see if a better fit is obtained by scaling old_skeleton up or down.
#' @return A named list containing elements new_skeleton (the updated skeleton) and new_scale
#' (the suggested revised prior scale to put on the single parameter beta in the power model).
#'
#' @references
#'
#' \insertRef{boonstra2020}{seamlesssim}
#'
#' \insertRef{dfcrm2019}{seamlesssim}
#'
#' @importFrom stats logLik glm binomial predict vcov
#' @export
calc_new_skeleton = function(old_skeleton,
                             tox,
                             level,
                             offset_seq = NULL) {

  if(is.null(offset_seq)) {
    offset_seq = exp(seq(-4, 4, length = 501));
  }
  log_likelihood_seq = rep(-Inf,length(offset_seq));
  fixed_covariate = log(-log(old_skeleton[level]));

  for(k in 1:length(offset_seq)) {
    log_likelihood_seq[k] = logLik(glm(I(1-tox) ~ offset(offset_seq[k]*fixed_covariate), family = binomial("cloglog")))
    if((k>2)&&(log_likelihood_seq[k]<log_likelihood_seq[k-1])) {break;}
  }
  best_offset = offset_seq[which.max(log_likelihood_seq)];
  best_model = glm(I(1 - tox) ~ offset(best_offset * fixed_covariate), family = binomial("cloglog"));
  list(new_skeleton = 1 - predict(best_model,newdata=data.frame(best_offset = best_offset, fixed_covariate = log(-log(old_skeleton))),type="response"),
       new_scale = as.numeric(sqrt(vcov(best_model))));
}
