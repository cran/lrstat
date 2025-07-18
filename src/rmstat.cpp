#include "utilities.h"

using namespace Rcpp;


//' @title Restricted Mean Survival Time
//' @description Obtains the restricted mean survival time over an interval.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//'
//' @return The integral of the survival function from \code{t1} to \code{t2}
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' rmst(t1 = 0, t2 = 7, piecewiseSurvivalTime = c(0, 6),
//'      lambda = c(0.0533, 0.0309))
//'
//' @export
// [[Rcpp::export]]
double rmst(const double t1 = 0,
            const double t2 = NA_REAL,
            const NumericVector& piecewiseSurvivalTime = 0,
            const NumericVector& lambda = NA_REAL) {

  if (std::isnan(t2)) {
    stop("t2 must be provided");
  }

  if (t1 < 0) {
    stop("t1 must be non-negative");
  }

  if (t2 < t1) {
    stop("t1 must be less than or equal to t2");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (piecewiseSurvivalTime.size() > 1 &&
      is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(lambda < 0))) {
    stop("lambda must be non-negative");
  }

  NumericVector t = piecewiseSurvivalTime;

  // identify the time interval containing the specified analysis time
  NumericVector time = NumericVector::create(t1, t2);
  IntegerVector m = findInterval3(time, t) - 1;

  double s1, s2, ch = 0.0, aval;
  for (int j=0; j<=m[0]-1; j++) {
    ch += lambda[j]*(t[j+1] - t[j]);
  }

  if (m[0] == m[1]) {
    s1 = exp(-lambda[m[0]]*(t1 - t[m[0]]));
    s2 = exp(-lambda[m[0]]*(t2 - t[m[0]]));
    aval = exp(-ch)*(s1 - s2)/lambda[m[0]];
  } else {
    s1 = exp(-lambda[m[0]]*(t1 - t[m[0]]));
    s2 = exp(-lambda[m[0]]*(t[m[0]+1] - t[m[0]]));
    aval = exp(-ch)*(s1 - s2)/lambda[m[0]];

    for (int j=m[0]+1; j<=m[1]-1; j++) {
      ch += lambda[j-1]*(t[j] - t[j-1]);
      s2 = exp(-lambda[j]*(t[j+1] - t[j]));
      aval = aval + exp(-ch)*(1.0 - s2)/lambda[j];
    }

    ch += lambda[m[1]-1]*(t[m[1]] - t[m[1]-1]);
    s2 = exp(-lambda[m[1]]*(t2 - t[m[1]]));
    aval = aval + exp(-ch)*(1.0 - s2)/lambda[m[1]];
  }

  return aval;
}


// define the integrand function for covrmst
struct rmparams {
  double time;
  double tau1;
  double tau2;
  double phi;
  NumericVector accrualTime;
  NumericVector accrualIntensity;
  NumericVector piecewiseSurvivalTime;
  NumericVector lambda;
  NumericVector gamma;
  double accrualDuration;
};


void f_rm(double *x, int n, void *ex) {
  rmparams *param = (rmparams *) ex;
  NumericVector u0(n), a1(n), a2(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
    a1[i] = rmst(u0[i], param->tau1, param->piecewiseSurvivalTime,
                 param->lambda);
    a2[i] = rmst(u0[i], param->tau2, param->piecewiseSurvivalTime,
                 param->lambda);
  }

  IntegerVector j = findInterval3(u0, param->piecewiseSurvivalTime) - 1;
  NumericVector lambda = param->lambda[j];
  NumericVector p = patrisk(u0, param->piecewiseSurvivalTime, param->lambda,
                            param->gamma);
  u0 = param->time - u0;
  NumericVector N = accrual(u0, param->accrualTime, param->accrualIntensity,
                            param->accrualDuration);
  u0 = a1*a2*lambda/((param->phi)*N*p);
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}



//' @title Covariance Between Restricted Mean Survival Times
//' @description Obtains the covariance between restricted mean survival
//' times at two different time points.
//'
//' @param t2 The calendar time for analysis 2.
//' @param tau1 The milestone time for analysis 1.
//' @param tau2 The milestone time for analysis 2.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_maxFollowupTime
//'
//' @return The covariance between the restricted mean survival times
//' for each treatment group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' covrmst(t2 = 25, tau1 = 16, tau2 = 18, allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericVector covrmst(const double t2 = NA_REAL,
                      const double tau1 = NA_REAL,
                      const double tau2 = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {

  if (std::isnan(t2)) {
    stop("t2 must be provided");
  }

  if (std::isnan(tau1)) {
    stop("tau1 must be provided");
  }

  if (tau1 <= 0) {
    stop("tau1 must be positive");
  }

  if (std::isnan(tau2)) {
    stop("tau2 must be provided");
  }

  if (tau2 < tau1) {
    stop("tau2 must be greater than or equal to tau1");
  }

  if (t2 <= tau2) {
    stop("t2 must be greater than tau2");
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (piecewiseSurvivalTime.size() > 1 &&
      is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(maxFollowupTime)) {
    stop("maxFollowupTime must be provided");
  }

  if (tau2 > maxFollowupTime) {
    stop("tau2 must be less than or equal to maxFollowupTime");
  }


  NumericVector t = piecewiseSurvivalTime;
  double phi = allocationRatioPlanned/(1 + allocationRatioPlanned);

  double tol = 1e-6;
  double umax = std::min(tau1, maxFollowupTime);

  rmparams param1 = {t2, tau1, tau2, phi, accrualTime, accrualIntensity,
                     t, lambda1, gamma1, accrualDuration};
  rmparams param2 = {t2, tau1, tau2, 1-phi, accrualTime, accrualIntensity,
                     t, lambda2, gamma2, accrualDuration};

  double q1 = quad(f_rm, &param1, 0.0, umax, tol)[0];
  double q2 = quad(f_rm, &param2, 0.0, umax, tol)[0];

  return NumericVector::create(q1, q2);
}


//' @title Restricted Mean Survival Time by Stratum
//' @description Obtains the restricted mean survival time and associated
//' variance by treatment group and by stratum at a given calendar time.
//'
//' @param time The calendar time for data cut.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//'
//' @return A data frame containing the following variables:
//'
//' * \code{stratum}: The stratum.
//'
//' * \code{time}: The calendar time since trial start.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{nevents}: The total number of events.
//'
//' * \code{nevents1}: The number of events in the active treatment group.
//'
//' * \code{nevents2}: The number of events in the control group.
//'
//' * \code{ndropouts}: The total number of dropouts.
//'
//' * \code{ndropouts1}: The number of dropouts in the active treatment
//'   group.
//'
//' * \code{ndropouts2}: The number of dropouts in the control group.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{nmilestone}: The total number of subjects reaching milestone.
//'
//' * \code{nmilestone1}: The number of subjects reaching milestone
//'   in the active treatment group.
//'
//' * \code{nmiletone2}: The number of subjects reaching milestone
//'   in the control group.
//'
//' * \code{rmst1}: The restricted mean survival time for the treatment
//'   group.
//'
//' * \code{rmst2}: The restricted mean survival time for the control group.
//'
//' * \code{rmstDiff}: The difference in restricted mean survival times,
//'   i.e., \code{rmst1 - rmst2}.
//'
//' * \code{vrmst1}: The variance for \code{rmst1}.
//'
//' * \code{vrmst2}: The variance for \code{rmst2}.
//'
//' * \code{vrmstDiff}: The variance for \code{rmstDiff}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' rmstat1(time = 40,
//'         milestone = 18,
//'         allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
//'         piecewiseSurvivalTime = c(0, 6),
//'         stratumFraction = c(0.2, 0.8),
//'         lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'         lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
DataFrame rmstat1(const double time = NA_REAL,
                  const double milestone = NA_REAL,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& lambda1 = NA_REAL,
                  const NumericVector& lambda2 = NA_REAL,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  const double accrualDuration = NA_REAL,
                  const double followupTime = NA_REAL,
                  const bool fixedFollowup = 0) {

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nintervals) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nintervals) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }


  // obtain the follow-up time for the first enrolled subject
  double maxFollowupTime;
  if (fixedFollowup) {
    maxFollowupTime = followupTime;
  } else {
    maxFollowupTime = accrualDuration + followupTime;
  }

  NumericVector ss(1, time);
  double a = accrual(ss, accrualTime, accrualIntensity, accrualDuration)[0];
  double phi = allocationRatioPlanned/(1 + allocationRatioPlanned);

  IntegerVector l1 = Range(0, nintervals-1);
  IntegerVector l(nintervals);
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericMatrix x(1,2), y(1,2);
  NumericVector nsubjects(nstrata);
  NumericMatrix nevents(nstrata, 2), ndropouts(nstrata, 2);
  NumericVector rmst1(nstrata), rmst2(nstrata), rmstDiff(nstrata);
  NumericVector vrmst1(nstrata), vrmst2(nstrata), vrmstDiff(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector calTime(nstrata), mileTime(nstrata);
  DataFrame df;

  NumericVector tt = NumericVector::create(time - milestone);
  double a2 = accrual(tt, accrualTime, accrualIntensity, accrualDuration)[0];
  NumericVector nmilestone1(nstrata), nmilestone2(nstrata);
  NumericVector nmilestone(nstrata);

  for (int h=0; h<nstrata; h++) {
    stratum[h] = h+1;
    calTime[h] = time;
    mileTime[h] = milestone;

    double frac = stratumFraction[h];
    l = h*nintervals + l1;
    lam1 = lambda1x[l];
    lam2 = lambda2x[l];
    gam1 = gamma1x[l];
    gam2 = gamma2x[l];

    // number of events in the stratum at the specified calendar time
    x = nevent2(ss, allocationRatioPlanned, accrualTime,
                frac*accrualIntensity,
                piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                accrualDuration, followupTime, maxFollowupTime);
    // number of dropouts in the stratum at the specified calendar time
    y = nevent2(ss, allocationRatioPlanned, accrualTime,
                frac*accrualIntensity,
                piecewiseSurvivalTime, gam1, gam2, lam1, lam2,
                accrualDuration, followupTime, maxFollowupTime);

    // obtain number of enrolled subjects and subjects having an event/dropout
    nsubjects[h] = frac*a;
    nevents(h, _) = x.row(0);
    ndropouts(h, _) = y.row(0);

    // obtain number of subjects reaching milestone
    double ncom = frac*a2;
    double p1 = patrisk(milestone, piecewiseSurvivalTime, lam1, gam1)[0];
    double p2 = patrisk(milestone, piecewiseSurvivalTime, lam2, gam2)[0];
    nmilestone1[h] = phi*ncom*p1;
    nmilestone2[h] = (1-phi)*ncom*p2;
    nmilestone[h] = nmilestone1[h] + nmilestone2[h];

    rmst1[h] = rmst(0, milestone, piecewiseSurvivalTime, lam1);
    rmst2[h] = rmst(0, milestone, piecewiseSurvivalTime, lam2);
    rmstDiff[h] = rmst1[h] - rmst2[h];

    NumericVector v = covrmst(
      time, milestone, milestone, allocationRatioPlanned,
      accrualTime, frac*accrualIntensity, piecewiseSurvivalTime,
      lam1, lam2, gam1, gam2, accrualDuration, maxFollowupTime);

    vrmst1[h] = v[0];
    vrmst2[h] = v[1];
    vrmstDiff[h] = v[0] + v[1];
  }

  // number of subjects having an event in each treatment group and overall
  NumericVector nevents1 = nevents(_, 0);
  NumericVector nevents2 = nevents(_, 1);
  NumericVector neventst = nevents1 + nevents2;

  NumericVector ndropouts1 = ndropouts(_, 0);
  NumericVector ndropouts2 = ndropouts(_, 1);
  NumericVector ndropoutst = ndropouts1 + ndropouts2;

  // output the requested information
  df = DataFrame::create(_["stratum"] = stratum,
                         _["time"] = calTime,
                         _["subjects"] = nsubjects,
                         _["nevents"] = neventst,
                         _["nevents1"] = nevents1,
                         _["nevents2"] = nevents2,
                         _["ndropouts"] = ndropoutst,
                         _["ndropouts1"] = ndropouts1,
                         _["ndropouts2"] = ndropouts2,
                         _["milestone"] = mileTime,
                         _["nmilestone"] = nmilestone,
                         _["nmilestone1"] = nmilestone1,
                         _["nmilestone2"] = nmilestone2,
                         _["rmst1"] = rmst1,
                         _["rmst2"] = rmst2,
                         _["rmstDiff"] = rmstDiff,
                         _["vrmst1"] = vrmst1,
                         _["vrmst2"] = vrmst2,
                         _["vrmstDiff"] = vrmstDiff);

  return df;
}


//' @title Stratified Difference in Restricted Mean Survival Times
//' @description Obtains the stratified restricted mean survival times
//' and difference in restricted mean survival times at given calendar
//' times.
//'
//' @param time A vector of calendar times for data cut.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//'
//' @return A data frame containing the following variables:
//'
//' * \code{time}: The calendar time since trial start.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{nevents}: The total number of events.
//'
//' * \code{nevents1}: The number of events in the active treatment group.
//'
//' * \code{nevents2}: The number of events in the control group.
//'
//' * \code{ndropouts}: The total number of dropouts.
//'
//' * \code{ndropouts1}: The number of dropouts in the active treatment
//'   group.
//'
//' * \code{ndropouts2}: The number of dropouts in the control group.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{nmilestone}: The total number of subjects reaching milestone.
//'
//' * \code{nmilestone1}: The number of subjects reaching milestone
//'   in the active treatment group.
//'
//' * \code{nmiletone2}: The number of subjects reaching milestone
//'   in the control group.
//'
//' * \code{rmst1}: The restricted mean survival time for the treatment
//'   group.
//'
//' * \code{rmst2}: The restricted mean survival time for the control group.
//'
//' * \code{rmstDiff}: The difference in restricted mean survival times,
//'   i.e., \code{rmst1 - rmst2}.
//'
//' * \code{vrmst1}: The variance for \code{rmst1}.
//'
//' * \code{vrmst2}: The variance for \code{rmst2}.
//'
//' * \code{vrmstDiff}: The variance for \code{rmstDiff}.
//'
//' * \code{information}: The information for \code{rmstDiff}, equal to
//'   \code{1/vrmstDiff}.
//'
//' * \code{rmstDiffZ}: The Z-statistic value, i.e.,
//'   \code{rmstDiff/sqrt(vrmstDiff)}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' rmstat(time = c(22, 40),
//'        milestone = 18,
//'        allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 8),
//'        accrualIntensity = 26/9*seq(1, 9),
//'        piecewiseSurvivalTime = c(0, 6),
//'        stratumFraction = c(0.2, 0.8),
//'        lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'        lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'        gamma1 = -log(1-0.05)/12,
//'        gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 22,
//'        followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
DataFrame rmstat(const NumericVector& time = NA_REAL,
                 const double milestone = NA_REAL,
                 const double allocationRatioPlanned = 1,
                 const NumericVector& accrualTime = 0,
                 const NumericVector& accrualIntensity = NA_REAL,
                 const NumericVector& piecewiseSurvivalTime = 0,
                 const NumericVector& stratumFraction = 1,
                 const NumericVector& lambda1 = NA_REAL,
                 const NumericVector& lambda2 = NA_REAL,
                 const NumericVector& gamma1 = 0,
                 const NumericVector& gamma2 = 0,
                 const double accrualDuration = NA_REAL,
                 const double followupTime = NA_REAL,
                 const bool fixedFollowup = 0) {

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);

  if (is_true(any(is_na(time)))) {
    stop("time must be provided");
  }

  if (is_true(any(time <= 0))) {
    stop("time must be positive");
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nintervals) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nintervals) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
    stop("followupTime must be provided");
  }

  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }

  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }

  if (is_true(any(time > accrualDuration + followupTime))) {
    stop("time cannot exceed accrualDuration + followupTime");
  }

  if (is_true(any(time <= milestone))) {
    stop("time must be greater than milestone");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }


  int k = static_cast<int>(time.size());
  NumericVector calTime(k), mileTime(k), subjects(k),
  nevents(k), nevents1(k), nevents2(k),
  ndropouts(k), ndropouts1(k), ndropouts2(k),
  nmilestone(k), nmilestone1(k), nmilestone2(k),
  rmst1(k), rmst2(k), vrmst1(k), vrmst2(k),
  rmstDiff(k), vrmstDiff(k), information(k), rmstDiffZ(k);
  DataFrame df;

  for (int j=0; j<k; j++) {
    df = rmstat1(time[j], milestone, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1x, lambda2x, gamma1x, gamma2x,
                 accrualDuration, followupTime, fixedFollowup);

    calTime[j] = max(NumericVector(df[1]));
    subjects[j] = sum(NumericVector(df[2]));
    nevents[j] = sum(NumericVector(df[3]));
    nevents1[j] = sum(NumericVector(df[4]));
    nevents2[j] = sum(NumericVector(df[5]));
    ndropouts[j] = sum(NumericVector(df[6]));
    ndropouts1[j] = sum(NumericVector(df[7]));
    ndropouts2[j] = sum(NumericVector(df[8]));
    mileTime[j] = max(NumericVector(df[9]));
    nmilestone[j] = sum(NumericVector(df[10]));
    nmilestone1[j] = sum(NumericVector(df[11]));
    nmilestone2[j] = sum(NumericVector(df[12]));
    rmst1[j] = sum(stratumFraction*NumericVector(df[13]));
    rmst2[j] = sum(stratumFraction*NumericVector(df[14]));
    rmstDiff[j] = rmst1[j] - rmst2[j];
    vrmst1[j] = sum(stratumFraction*stratumFraction*NumericVector(df[16]));
    vrmst2[j] = sum(stratumFraction*stratumFraction*NumericVector(df[17]));
    vrmstDiff[j] = vrmst1[j] + vrmst2[j];
    information[j] = 1.0/vrmstDiff[j];
    rmstDiffZ[j] = rmstDiff[j]*sqrt(information[j]);
  }

  df = DataFrame::create(_["time"] = calTime,
                         _["subjects"] = subjects,
                         _["nevents"] = nevents,
                         _["nevents1"] = nevents1,
                         _["nevents2"] = nevents2,
                         _["ndropouts"] = ndropouts,
                         _["ndropouts1"] = ndropouts1,
                         _["ndropouts2"] = ndropouts2,
                         _["milestone"] = mileTime,
                         _["nmilestone"] = nmilestone,
                         _["nmilestone1"] = nmilestone1,
                         _["nmilestone2"] = nmilestone2,
                         _["rmst1"] = rmst1,
                         _["rmst2"] = rmst2,
                         _["rmstDiff"] = rmstDiff,
                         _["vrmst1"] = vrmst1,
                         _["vrmst2"] = vrmst2,
                         _["vrmstDiff"] = vrmstDiff,
                         _["information"] = information,
                         _["rmstDiffZ"] = rmstDiffZ);

  return df;
}


//' @title Power for Difference in Restricted Mean Survival Times
//' @description Estimates the power for testing the difference in
//' restricted mean survival times in a two-sample survival design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{rmpower} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{rmstDiffH0}: The difference in restricted mean survival
//'       times under the null hypothesis.
//'
//'     - \code{rmst1}: The restricted mean survival time for the
//'       treatment group.
//'
//'     - \code{rmst2}: The restricted mean survival time for the
//'       control group.
//'
//'     - \code{rmstDiff}: The difference in restricted mean survival times,
//'       equal to \code{rmst1 - rmst2}.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRmstDiff}: The efficacy boundaries on the restricted
//'       mean survival time difference scale.
//'
//'     - \code{futilityRmstDiff}: The futility boundaries on the restricted
//'       mean survival time difference scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   and \code{spendingTime}.
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by
//'   treatment group:
//'
//'     - \code{numberOfEvents1}: The number of events by stage for
//'       the treatment group.
//'
//'     - \code{numberOfDropouts1}: The number of dropouts by stage for
//'       the treatment group.
//'
//'     - \code{numberOfSubjects1}: The number of subjects by stage for
//'       the treatment group.
//'
//'     - \code{numberOfMilestone1}: The number of subjects reaching
//'       milestone by stage for the active treatment group.
//'
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{numberOfMilestone2}: The number of subjects reaching
//'       milestone by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfMilestone1}: The expected number of subjects
//'       reaching milestone for the active treatment group.
//'
//'     - \code{expectedNumberOfEvents2}: The expected number of events for
//'       control group.
//'
//'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
//'       for the control group.
//'
//'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
//'       for the control group.
//'
//'     - \code{expectedNumberOfMilestone2}: The expected number of subjects
//'       reaching milestone for the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' rmpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         milestone = 18,
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'         accrualIntensity = 100/9*seq(1, 9),
//'         piecewiseSurvivalTime = c(0, 6),
//'         stratumFraction = c(0.2, 0.8),
//'         lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'         lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List rmpower(const int kMax = 1,
             const NumericVector& informationRates = NA_REAL,
             const LogicalVector& efficacyStopping = NA_LOGICAL,
             const LogicalVector& futilityStopping = NA_LOGICAL,
             const NumericVector& criticalValues = NA_REAL,
             const double alpha = 0.025,
             const std::string typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& userAlphaSpending = NA_REAL,
             const NumericVector& futilityBounds = NA_REAL,
             const std::string typeBetaSpending = "none",
             const double parameterBetaSpending = NA_REAL,
             const double milestone = NA_REAL,
             const double rmstDiffH0 = 0,
             const double allocationRatioPlanned = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = NA_REAL,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const NumericVector& lambda1 = NA_REAL,
             const NumericVector& lambda2 = NA_REAL,
             const NumericVector& gamma1 = 0,
             const NumericVector& gamma2 = 0,
             const double accrualDuration = NA_REAL,
             const double followupTime = NA_REAL,
             const bool fixedFollowup = 0,
             const NumericVector& spendingTime = NA_REAL,
             const double studyDuration = NA_REAL) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }

    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop("parameterBetaSpending must be positive for sfKD");
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (lambda1.size() != 1 && lambda1.size() != nintervals &&
      lambda1.size() != nsi) {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() != 1 && lambda2.size() != nintervals &&
      lambda2.size() != nsi) {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
    stop("followupTime must be provided");
  }

  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }

  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));
  alpha1 = cumAlphaSpent[kMax - 1];

  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds && futilityBounds.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    } else if (!missingFutilityBounds && futilityBounds.size() < kMax-1) {
      stop("Insufficient length of futilityBounds");
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }


  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  // obtain the timing of interim analysis
  DataFrame rm;
  NumericVector time(kMax);
  NumericVector u0(1, studyDuration1);
  rm = rmstat(u0, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double maxInformation = sum(NumericVector(rm[18]));
  double rmst1 = sum(NumericVector(rm[12]));
  double rmst2 = sum(NumericVector(rm[13]));
  double rmstDiff = sum(NumericVector(rm[14]));
  NumericVector theta(kMax, rmstDiff);

  // information at milestone
  u0[0] = milestone + 1.0e-6;
  rm = rmstat(u0, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double information1 = sum(NumericVector(rm[18]));
  if (informationRates1[0] <= information1/maxInformation) {
    std::string str1 = "The first information rate must exceed that";
    std::string str2 = "at the milestone time\n";
    std::string str3 = "The information at the milestone time is";
    std::string str4 = "The maximum information at study end is";
    std::string errmsg = str1 + " " + str2 +
      str3 + " " + std::to_string(information1) + "\n" +
      str4 + " " + std::to_string(maxInformation) + "\n";
    stop(errmsg);
  }

  NumericVector I = maxInformation*informationRates1;


  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup);
              return sum(NumericVector(rm[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    information1 = I[i];
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  };
  time[kMax-1] = studyDuration1;

  rm = rmstat(time, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double phi = allocationRatioPlanned/(allocationRatioPlanned+1);
  NumericVector nsubjects = NumericVector(rm[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
  NumericVector nevents = NumericVector(rm[2]);
  NumericVector nevents1 = NumericVector(rm[3]);
  NumericVector nevents2 = NumericVector(rm[4]);
  NumericVector ndropouts = NumericVector(rm[5]);
  NumericVector ndropouts1 = NumericVector(rm[6]);
  NumericVector ndropouts2 = NumericVector(rm[7]);
  NumericVector nmilestone = NumericVector(rm[9]);
  NumericVector nmilestone1 = NumericVector(rm[10]);
  NumericVector nmilestone2 = NumericVector(rm[11]);

  // compute the stagewise exit probabilities for efficacy and futility
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprobcpp(criticalValues1, futilityBounds1, theta, I);
  } else {
    NumericVector w(kMax, 1.0);
    List out = getPower(alpha1, kMax, criticalValues1, theta, I,
                        bsf, bsfpar, spendingTime1, futilityStopping1, w);
    futilityBounds1 = out[1];
    probs = out[2];
  }

  NumericVector efficacyP(kMax);
  NumericVector futilityP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
    futilityP[i] = 1 - R::pnorm(futilityBounds1[i], 0, 1, 1, 0);
  }

  // stagewise total exit probabilities
  NumericVector pu = NumericVector(probs[0]);
  NumericVector pl = NumericVector(probs[1]);
  NumericVector ptotal = pu + pl;

  double overallReject = sum(pu);
  double expectedNumberOfEvents = sum(ptotal*nevents);
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedNumberOfEvents1 = sum(ptotal*nevents1);
  double expectedNumberOfDropouts1 = sum(ptotal*ndropouts1);
  double expectedNumberOfSubjects1 = sum(ptotal*nsubjects1);
  double expectedNumberOfMilestone1 = sum(ptotal*nmilestone1);
  double expectedNumberOfEvents2 = sum(ptotal*nevents2);
  double expectedNumberOfDropouts2 = sum(ptotal*ndropouts2);
  double expectedNumberOfSubjects2 = sum(ptotal*nsubjects2);
  double expectedNumberOfMilestone2 = sum(ptotal*nmilestone2);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);

  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  NumericVector rdu = rmstDiffH0 + criticalValues1/sqrt(I);
  NumericVector rdl = rmstDiffH0 + futilityBounds1/sqrt(I);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      rdu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      rdl[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }


  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["numberOfEvents"] = nevents,
    _["numberOfDropouts"] = ndropouts,
    _["numberOfSubjects"] = nsubjects,
    _["numberOfMilestone"] = nmilestone,
    _["analysisTime"] = time,
    _["efficacyRmstDiff"] = rdu,
    _["futilityRmstDiff"] = rdl,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["rmstDiffH0"] = rmstDiffH0,
    _["rmst1"] = rmst1,
    _["rmst2"] = rmst2,
    _["rmstDiff"] = rmstDiff);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["stratumFraction"] = stratumFraction,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["spendingTime"] = spendingTime);

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["numberOfMilestone1"] = nmilestone1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["numberOfMilestone2"] = nmilestone2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedNumberOfMilestone1"] = expectedNumberOfMilestone1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2,
    _["expectedNumberOfMilestone2"] = expectedNumberOfMilestone2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "rmpower";

  return result;
}


//' @title Sample Size for Difference in Restricted Mean Survival Times
//' @description Obtains the needed accrual duration given power,
//' accrual intensity, and follow-up time, the needed follow-up time
//' given power, accrual intensity, and accrual duration, or the needed
//' absolute accrual intensity given power, relative accrual intensity,
//' accrual duration, and follow-up time in a two-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupTime, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{rmpower} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{rmpower} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpower}}
//'
//' @examples
//' # Example 1: Obtains follow-up time given power, accrual intensity,
//' # and accrual duration for variable follow-up. Of note, the power
//' # reaches the maximum when the follow-up time equals milestone.
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up time for variable follow-up
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up time for fixed follow-up
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = NA,
//'              followupTime = 18, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List rmsamplesize(const double beta = 0.2,
                  const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const LogicalVector& efficacyStopping = NA_LOGICAL,
                  const LogicalVector& futilityStopping = NA_LOGICAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.025,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const NumericVector& futilityBounds = NA_REAL,
                  const std::string typeBetaSpending = "none",
                  const double parameterBetaSpending = NA_REAL,
                  const NumericVector& userBetaSpending = NA_REAL,
                  const double milestone = NA_REAL,
                  const double rmstDiffH0 = 0,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& lambda1 = NA_REAL,
                  const NumericVector& lambda2 = NA_REAL,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  double accrualDuration = NA_REAL,
                  double followupTime = NA_REAL,
                  const bool fixedFollowup = 0,
                  const NumericVector& interval =
                    NumericVector::create(0.001, 240),
                    const NumericVector& spendingTime = NA_REAL,
                    const bool rounding = 1) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi);


  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }

    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(futilityBounds))) && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }

  if (!std::isnan(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }

  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }

    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }

  if (fixedFollowup && std::isnan(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (!std::isnan(accrualDuration) && !std::isnan(followupTime) &&
      (milestone >= accrualDuration + followupTime)) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (interval.size() != 2) {
    stop("interval must have 2 elements");
  }

  if (interval[0] < 0) {
    stop("lower limit of interval must be positive");
  }

  if (interval[0] >= interval[1]) {
    stop("upper limit must be greater than lower limit for interval");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds && futilityBounds.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    } else if (!missingFutilityBounds && futilityBounds.size() < kMax-1) {
      stop("Insufficient length of futilityBounds");
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }


  std::string unknown;
  // search for the solution according to the input
  if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualDuration";
  } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
    unknown = "followupTime";
  } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }


  NumericVector rmsts1(nstrata), rmsts2(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  for (int h=0; h<nstrata; h++) {
    l = h*nintervals + l1;
    NumericVector lam1 = lambda1x[l];
    NumericVector lam2 = lambda2x[l];
    rmsts1[h] = rmst(0, milestone, piecewiseSurvivalTime, lam1);
    rmsts2[h] = rmst(0, milestone, piecewiseSurvivalTime, lam2);
  }

  double rmst1 = sum(stratumFraction*rmsts1);
  double rmst2 = sum(stratumFraction*rmsts2);
  double theta1 = rmst1 - rmst2 - rmstDiffH0;
  NumericVector theta(kMax, theta1);

  List design = getDesign(
    beta, NA_REAL, theta1, kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
    bsf, bsfpar, userBetaSpending, spendingTime1, 1);

  DataFrame byStageResults = DataFrame(design["byStageResults"]);
  futilityBounds1 = byStageResults["futilityBounds"];

  DataFrame overallResults = DataFrame(design["overallResults"]);
  double maxInformation = overallResults["information"];
  double studyDuration;

  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            unknown, maxInformation](double aval)-> double{
              NumericVector accrualIntensity1 = clone(accrualIntensity);
              double dur1=0, dur2=0;

              if (unknown == "accrualDuration") {
                dur1 = aval;
                dur2 = followupTime;
              } else if (unknown == "followupTime") {
                dur1 = accrualDuration;
                dur2 = aval;
              } else if (unknown == "accrualIntensity") {
                dur1 = accrualDuration;
                dur2 = followupTime;
                accrualIntensity1 = aval*accrualIntensity;
              }

              // obtain the maximum information at study end
              NumericVector u0(1, dur1 + dur2);
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup);
              return sum(NumericVector(rm[18])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
    accrualDuration = brent(f, lower, interval[1], 1.0e-6);
    studyDuration = accrualDuration + followupTime;
  } else if (unknown == "followupTime") {
    if (f(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[18])) - maxInformation;
                };

      accrualDuration = brent(g, accrualDuration, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // adjust follow-up time
      double lower =  std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(f, lower, milestone, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  } else if (unknown == "accrualIntensity") {
    if (accrualDuration + followupTime <= milestone) {
      std::string str1 = "NOTE: The followupTime is too short for ";
      std::string str2 = "anyone to reach milestone.";
      std::string str3 = "NOTE: followupTime is increased to ";
      std::string str4 = "milestone.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";
      followupTime = milestone;
    }
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
    studyDuration = accrualDuration + followupTime;
  }


  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup);
                    return sum(NumericVector(rm[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // the follow-up time cannot be adjusted for fixed follow-up
        // adjust study duration to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup);
                    return sum(NumericVector(rm[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }

  resultH1 = rmpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  overallResults = DataFrame(resultH1["overallResults"]);
  maxInformation = overallResults["information"];

  // obtain results under H0 by matching the maximum information
  // first find the hazard rate for the active treatment group that yields
  // the specified difference in restricted mean survival times under H0
  auto frmst = [milestone, piecewiseSurvivalTime, stratumFraction,
                nintervals, nstrata, l1, lambda2x, rmst2,
                rmstDiffH0](double aval)-> double {
                  NumericVector rmsts1(nstrata);
                  for (int h=0; h<nstrata; h++) {
                    IntegerVector l = h*nintervals + l1;
                    NumericVector lam2 = lambda2x[l];
                    rmsts1[h] = rmst(0, milestone, piecewiseSurvivalTime,
                                     aval*lam2);
                  }
                  double rmst1 = sum(stratumFraction*rmsts1);
                  return rmst1 - rmst2 - rmstDiffH0;
                };

  double lower = 0.001, upper = 1.999;
  while (frmst(upper) > 0) {
    lower = upper;
    upper = 2.0*upper;
  }
  double aval = brent(frmst, lower, upper, 1.0e-6);
  NumericVector lambda1H0 = aval*lambda2;

  if (!fixedFollowup) {
    auto h = [milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda1H0, lambda2, gamma1, gamma2,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame rm = rmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  accrualDuration, aval, fixedFollowup);
                return sum(NumericVector(rm[18])) - maxInformation;
              };

    if (h(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1H0, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[18])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
    } else if (accrualDuration <= milestone || h(0) < 0) {
      // adjust follow-up time
      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(h, lower, milestone, 1.0e-6);
    } else {
      // decrease accrual duration
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1H0, lambda2, gamma1, gamma2,
                fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    aval, 0, fixedFollowup);
                  return sum(NumericVector(rm[18])) - maxInformation;
                };

      double lower = milestone + 1.0e-6;
      accrualDuration = brent(g, lower, accrualDuration, 1.0e-6);
      followupTime = 0.0;
    }
    studyDuration = accrualDuration + followupTime;
  } else { // fixed follow-up
    // cannot adjust the followupTime
    auto h = [milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda1H0, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame rm = rmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup);
                return sum(NumericVector(rm[18])) - maxInformation;
              };

    auto g = [milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda1H0, lambda2, gamma1, gamma2,
              followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, aval + followupTime);
                DataFrame rm = rmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  aval, followupTime, fixedFollowup);
                return sum(NumericVector(rm[18])) - maxInformation;
              };

    if (h(followupTime) < 0) { // increase accrual duration
      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (accrualDuration <= milestone || h(0) < 0) {
      // decrease study duration
      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double upper = std::min(milestone, followupTime);
      double aval = brent(h, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + aval;
    } else { // decrease accrual duration
      double lower = milestone + 1.0e-6;
      accrualDuration = brent(g, lower, accrualDuration, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = rmpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1H0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}



//' @title Power for One-Sample Restricted Mean Survival Time
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a one-group survival design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstH0 The restricted mean survival time under the null
//'   hypothesis.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda A vector of hazard rates for the event in each analysis
//'  time interval by stratum under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{rmpower1s} object with 3 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{numberOfMilestone}: The total number of subjects reaching
//'       milestone.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedNumberOfMilestone}: The expected number of subjects
//'       reaching milestone.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time to calculate the restricted
//'       mean survival time.
//'
//'     - \code{rmstH0}: The restricted mean survival time under the null
//'       hypothesis.
//'
//'     - \code{rmst}: The restricted mean survival time under the
//'       alternative hypothesis.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRmst}: The efficacy boundaries on the restricted
//'       mean survival time.
//'
//'     - \code{futilityRmst}: The futility boundaries on the restricted
//'       mean survival time.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{accrualTime},
//'   \code{accuralIntensity}, \code{piecewiseSurvivalTime},
//'   \code{stratumFraction}, \code{lambda}, \code{gamma},
//'   and \code{spendingTime}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' rmpower1s(kMax = 2, informationRates = c(0.8, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           milestone = 18, rmstH0 = 10,
//'           accrualTime = seq(0, 8),
//'           accrualIntensity = 26/9*seq(1, 9),
//'           piecewiseSurvivalTime = c(0, 6),
//'           stratumFraction = c(0.2, 0.8),
//'           lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'           gamma = -log(1-0.05)/12, accrualDuration = 22,
//'           followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
//'
// [[Rcpp::export]]
List rmpower1s(const int kMax = 1,
               const NumericVector& informationRates = NA_REAL,
               const LogicalVector& efficacyStopping = NA_LOGICAL,
               const LogicalVector& futilityStopping = NA_LOGICAL,
               const NumericVector& criticalValues = NA_REAL,
               const double alpha = 0.025,
               const std::string typeAlphaSpending = "sfOF",
               const double parameterAlphaSpending = NA_REAL,
               const NumericVector& userAlphaSpending = NA_REAL,
               const NumericVector& futilityBounds = NA_REAL,
               const std::string typeBetaSpending = "none",
               const double parameterBetaSpending = NA_REAL,
               const double milestone = NA_REAL,
               const double rmstH0 = NA_REAL,
               const NumericVector& accrualTime = 0,
               const NumericVector& accrualIntensity = NA_REAL,
               const NumericVector& piecewiseSurvivalTime = 0,
               const NumericVector& stratumFraction = 1,
               const NumericVector& lambda = NA_REAL,
               const NumericVector& gamma = 0,
               const double accrualDuration = NA_REAL,
               const double followupTime = NA_REAL,
               const bool fixedFollowup = 0,
               const NumericVector& spendingTime = NA_REAL,
               const double studyDuration = NA_REAL) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }

    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(rmstH0)) {
    stop("rmstH0 must be provided");
  }

  if (rmstH0 <= 0) {
    stop("rmstH0 must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(lambda < 0))) {
    stop("lambda must be non-negative");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (lambda.size() != 1 && lambda.size() != nintervals &&
      lambda.size() != nsi) {
    stop("Invalid length for lambda");
  }

  if (gamma.size() != 1 && gamma.size() != nintervals) {
    stop("Invalid length for gamma");
  }

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
    stop("followupTime must be provided");
  }

  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }

  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));
  alpha1 = cumAlphaSpent[kMax - 1];

  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds && futilityBounds.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    } else if (!missingFutilityBounds && futilityBounds.size() < kMax-1) {
      stop("Insufficient length of futilityBounds");
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }


  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }


  // obtain the timing of interim analysis using the twin treatment group
  DataFrame rm;
  NumericVector time(kMax);
  NumericVector u0(1, studyDuration1);
  rm = rmstat(u0, milestone, 1,
              accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup);

  double maxInformation = 2.0*sum(NumericVector(rm[18]));
  double rmstH1 = sum(NumericVector(rm[12]));
  NumericVector theta = rep(rmstH1 - rmstH0, kMax);

  // information at milestone
  u0[0] = milestone + 1.0e-6;
  rm = rmstat(u0, milestone, 1,
              accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup);

  double information1 = sum(NumericVector(rm[18]));
  if (informationRates1[0] <= information1/maxInformation) {
    std::string str1 = "The first information rate must exceed that";
    std::string str2 = "at the milestone time\n";
    std::string str3 = "The information at the milestone time is";
    std::string str4 = "The maximum information at study end is";
    std::string errmsg = str1 + " " + str2 +
      str3 + " " + std::to_string(information1) + "\n" +
      str4 + " " + std::to_string(maxInformation) + "\n";
    stop(errmsg);
  }

  NumericVector I = maxInformation*informationRates1;


  auto f = [milestone, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame rm = rmstat(
                u0, milestone, 1,
                accrualTime, 2.0*accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda, lambda, gamma, gamma,
                accrualDuration, followupTime, fixedFollowup);
              return 2.0*sum(NumericVector(rm[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information1 = std::max(I[i], 0.0);
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  }
  time[kMax-1] = studyDuration1;

  rm = rmstat(time, milestone, 1,
              accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup);

  NumericVector nsubjects = 0.5*NumericVector(rm[1]);
  NumericVector nevents = NumericVector(rm[3]);
  NumericVector ndropouts = NumericVector(rm[6]);
  NumericVector nmilestone = NumericVector(rm[10]);

  // compute the stagewise exit probabilities for efficacy and futility
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprobcpp(criticalValues1, futilityBounds1, theta, I);
  } else {
    NumericVector w(kMax, 1.0);
    List out = getPower(alpha1, kMax, criticalValues1, theta, I,
                        bsf, bsfpar, spendingTime1, futilityStopping1, w);
    futilityBounds1 = out[1];
    probs = out[2];
  }

  NumericVector efficacyP(kMax);
  NumericVector futilityP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
    futilityP[i] = 1 - R::pnorm(futilityBounds1[i], 0, 1, 1, 0);
  }

  // stagewise total exit probabilities
  NumericVector pu = NumericVector(probs[0]);
  NumericVector pl = NumericVector(probs[1]);
  NumericVector ptotal = pu + pl;

  double overallReject = sum(pu);
  double expectedNumberOfEvents = sum(ptotal*nevents);
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedNumberOfMilestone = sum(ptotal*nmilestone);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);

  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  NumericVector rmstu = rmstH0 + criticalValues1/sqrt(I);
  NumericVector rmstl = rmstH0 + futilityBounds1/sqrt(I);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      rmstu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      rmstl[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }


  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["numberOfEvents"] = nevents,
    _["numberOfDropouts"] = ndropouts,
    _["numberOfSubjects"] = nsubjects,
    _["numberOfMilestone"] = nmilestone,
    _["analysisTime"] = time,
    _["efficacyRmst"] = rmstu,
    _["futilityRmst"] = rmstl,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["numberOfMilestone"] = (nmilestone[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedNumberOfMilestone"] = expectedNumberOfMilestone,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["rmstH0"] = rmstH0,
    _["rmst"] = rmstH1);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["stratumFraction"] = stratumFraction,
    _["lambda"] = lambda,
    _["gamma"] = gamma,
    _["spendingTime"] = spendingTime);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "rmpower1s";

  return result;
}


//' @title Sample Size for One-Sample Restricted Mean Survival Time
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up duration, and relative accrual
//' rates in a one-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted survival time.
//' @param rmstH0 The restricted mean survival time under the null
//'   hypothesis.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda A vector of hazard rates for the event in each analysis
//'  time interval by stratum under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{rmpower1s} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{rmpower1s} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpower1s}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = 22,
//'                followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up duration for variable follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = 22,
//'                followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up duration for fixed follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = NA,
//'                followupTime = 18, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List rmsamplesize1s(const double beta = 0.2,
                    const int kMax = 1,
                    const NumericVector& informationRates = NA_REAL,
                    const LogicalVector& efficacyStopping = NA_LOGICAL,
                    const LogicalVector& futilityStopping = NA_LOGICAL,
                    const NumericVector& criticalValues = NA_REAL,
                    const double alpha = 0.025,
                    const std::string typeAlphaSpending = "sfOF",
                    const double parameterAlphaSpending = NA_REAL,
                    const NumericVector& userAlphaSpending = NA_REAL,
                    const NumericVector& futilityBounds = NA_REAL,
                    const std::string typeBetaSpending = "none",
                    const double parameterBetaSpending = NA_REAL,
                    const NumericVector& userBetaSpending = NA_REAL,
                    const double milestone = NA_REAL,
                    const double rmstH0 = NA_REAL,
                    const NumericVector& accrualTime = 0,
                    const NumericVector& accrualIntensity = NA_REAL,
                    const NumericVector& piecewiseSurvivalTime = 0,
                    const NumericVector& stratumFraction = 1,
                    const NumericVector& lambda = NA_REAL,
                    const NumericVector& gamma = 0,
                    double accrualDuration = NA_REAL,
                    double followupTime = NA_REAL,
                    const bool fixedFollowup = 0,
                    const NumericVector& interval =
                      NumericVector::create(0.001, 240),
                      const NumericVector& spendingTime = NA_REAL,
                      const bool rounding = 1) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector lambdax(nsi);


  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }

    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(futilityBounds))) && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(rmstH0)) {
    stop("rmstH0 must be provided");
  }

  if (rmstH0 <= 0) {
    stop("rmstH0 must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(lambda < 0))) {
    stop("lambda must be non-negative");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (lambda.size() == 1) {
    lambdax = rep(lambda, nsi);
  } else if (lambda.size() == nintervals) {
    lambdax = rep(lambda, nstrata);
  } else if (lambda.size() == nsi) {
    lambdax = lambda;
  } else {
    stop("Invalid length for lambda");
  }

  if (lambda.size() != 1 && lambda.size() != nintervals &&
      lambda.size() != nsi) {
    stop("Invalid length for lambda");
  }

  if (gamma.size() != 1 && gamma.size() != nintervals) {
    stop("Invalid length for gamma");
  }

  if (!std::isnan(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }

  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }

    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }

  if (fixedFollowup && std::isnan(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (!std::isnan(accrualDuration) && !std::isnan(followupTime) &&
      (milestone >= accrualDuration + followupTime)) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (interval.size() != 2) {
    stop("interval must have 2 elements");
  }

  if (interval[0] < 0) {
    stop("lower limit of interval must be positive");
  }

  if (interval[0] >= interval[1]) {
    stop("upper limit must be greater than lower limit for interval");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds && futilityBounds.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    } else if (!missingFutilityBounds && futilityBounds.size() < kMax-1) {
      stop("Insufficient length of futilityBounds");
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }


  std::string unknown;
  // search for the solution according to the input
  if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualDuration";
  } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
    unknown = "followupTime";
  } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }


  NumericVector rmsts(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  for (int h=0; h<nstrata; h++) {
    l = h*nintervals + l1;
    NumericVector lam = lambdax[l];
    rmsts[h] = rmst(0, milestone, piecewiseSurvivalTime, lam);
  }

  double rmstH1 = sum(stratumFraction*rmsts);
  double theta1 = rmstH1 - rmstH0;
  NumericVector theta(kMax, theta1);

  List design = getDesign(
    beta, NA_REAL, theta1, kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
    bsf, bsfpar, userBetaSpending, spendingTime1, 1);

  DataFrame byStageResults = DataFrame(design["byStageResults"]);
  futilityBounds1 = byStageResults["futilityBounds"];

  DataFrame overallResults = DataFrame(design["overallResults"]);
  double maxInformation = overallResults["information"];
  double studyDuration;

  auto f = [milestone, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            unknown, maxInformation](double aval)-> double{
              NumericVector accrualIntensity1 = clone(accrualIntensity);
              double dur1=0, dur2=0;

              if (unknown == "accrualDuration") {
                dur1 = aval;
                dur2 = followupTime;
              } else if (unknown == "followupTime") {
                dur1 = accrualDuration;
                dur2 = aval;
              } else if (unknown == "accrualIntensity") {
                dur1 = accrualDuration;
                dur2 = followupTime;
                accrualIntensity1 = aval*accrualIntensity;
              }

              // obtain the maximum information at study end
              NumericVector u0(1, dur1 + dur2);
              DataFrame rm = rmstat(
                u0, milestone, 1,
                accrualTime, 2.0*accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda, lambda, gamma, gamma,
                dur1, dur2, fixedFollowup);
              return 2.0*sum(NumericVector(rm[18])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
    accrualDuration = brent(f, lower, interval[1], 1.0e-6);
    studyDuration = accrualDuration + followupTime;
  } else if (unknown == "followupTime") {
    if (f(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda, gamma,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda, lambda, gamma, gamma,
                    aval, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[18])) - maxInformation;
                };

      accrualDuration = brent(g, accrualDuration, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // adjust follow-up time
      double lower =  std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(f, lower, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  } else if (unknown == "accrualIntensity") {
    if (accrualDuration + followupTime <= milestone) {
      std::string str1 = "NOTE: The followupTime is too short for ";
      std::string str2 = "anyone to reach milestone.";
      std::string str3 = "NOTE: followupTime is increased to ";
      std::string str4 = "milestone.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";
      followupTime = milestone;
    }
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
    studyDuration = accrualDuration + followupTime;
  }


  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [milestone, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda, gamma,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, 1,
                      accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, aval, fixedFollowup);
                    return 2.0*sum(NumericVector(rm[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // the follow-up time cannot be adjusted for fixed follow-up
        // adjust study duration to obtain the target maximum information
        auto h = [milestone, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda, gamma,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, 1,
                      accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, followupTime, fixedFollowup);
                    return 2.0*sum(NumericVector(rm[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }

  resultH1 = rmpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  overallResults = DataFrame(resultH1["overallResults"]);
  maxInformation = overallResults["information"];


  // obtain results under H0 by matching the maximum information
  // first find the hazard rate that yields
  // the specified milestone survival probability under H0
  auto frmst = [milestone, piecewiseSurvivalTime, stratumFraction,
                nintervals, nstrata, l1, lambdax,
                rmstH0](double aval)-> double {
                  NumericVector rmsts(nstrata);
                  for (int h=0; h<nstrata; h++) {
                    IntegerVector l = h*nintervals + l1;
                    NumericVector lam = lambdax[l];
                    rmsts[h] = rmst(0, milestone, piecewiseSurvivalTime,
                                    aval*lam);
                  }
                  double rmstH1 = sum(stratumFraction*rmsts);
                  return rmstH1 - rmstH0;
                };

  double lower = 0.001, upper = 1.999;
  while (frmst(upper) > 0) {
    lower = upper;
    upper = 2.0*upper;
  }
  double aval = brent(frmst, lower, upper, 1.0e-6);
  NumericVector lambdaH0 = aval*lambda;

  if (!fixedFollowup) {
    auto h = [milestone, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambdaH0, gamma,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame rm = rmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  accrualDuration, aval, fixedFollowup);
                return 2.0*sum(NumericVector(rm[18])) - maxInformation;
              };

    if (h(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambdaH0, gamma,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[18])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
    } else if (accrualDuration <= milestone || h(0) < 0) {
      // adjust follow-up time
      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(h, lower, milestone, 1.0e-6);
    } else {
      // adjust accrual duration
      auto g = [milestone, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambdaH0, gamma, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval, 0, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[18])) - maxInformation;
                };

      double lower = milestone + 1.0e-6;
      accrualDuration = brent(g, lower, accrualDuration, 1.0e-6);
      followupTime = 0.0;
    }
    studyDuration = accrualDuration + followupTime;
  } else { // fixed follow-up
    // cannot adjust the followupTime
    auto h = [milestone, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambdaH0, gamma,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame rm = rmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  accrualDuration, followupTime, fixedFollowup);
                return 2.0*sum(NumericVector(rm[18])) - maxInformation;
              };

    auto g = [milestone, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambdaH0, gamma,
              followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, aval + followupTime);
                DataFrame rm = rmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  aval, followupTime, fixedFollowup);
                return 2.0*sum(NumericVector(rm[18])) - maxInformation;
              };

    if (h(followupTime) < 0) { // increase accrual duration
      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (accrualDuration <= milestone || h(0) < 0) {
      // decrease study duration
      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double upper = std::min(milestone, followupTime);
      double aval = brent(h, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + aval;
    } else { // decrease accrual duration
      double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
      accrualDuration = brent(g, lower, accrualDuration, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = rmpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambdaH0, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}



//' @title Power for Equivalence in Restricted Mean Survival Time Difference
//' @description Obtains the power for equivalence in restricted mean
//' survival time difference.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffLower The lower equivalence limit of restricted mean
//'   survival time difference.
//' @param rmstDiffUpper The upper equivalence limit of restricted mean
//'   survival time difference.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{rmpowerequiv} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{rmstDiffLower}: The lower equivalence limit of restricted
//'       mean survival time difference.
//'
//'     - \code{rmstDiffUpper}: The upper equivalence limit of restricted
//'       mean survival time difference.
//'
//'     - \code{rmst1}: The restricted mean survival time for the
//'       treatment group.
//'
//'     - \code{rmst2}: The restricted mean survival time for the
//'       control group.
//'
//'     - \code{rmstDiff}: The restricted mean survival time difference.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
//'       each of the two one-sided tests.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
//'       the two one-sided tests.
//'
//'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
//'       under \code{H10}.
//'
//'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
//'       under \code{H20}.
//'
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRmstDiffLower}: The efficacy boundaries on the
//'       restricted mean survival time difference scale for the one-sided
//'       null hypothesis at the lower equivalence limit.
//'
//'     - \code{efficacyRmstDiffUpper}: The efficacy boundaries on the
//'       restricted mean survival time difference scale for the one-sided
//'       null hypothesis at the upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   and \code{spendingTime}.
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by
//'   treatment group:
//'
//'     - \code{numberOfEvents1}: The number of events by stage for
//'       the treatment group.
//'
//'     - \code{numberOfDropouts1}: The number of dropouts by stage for
//'       the treatment group.
//'
//'     - \code{numberOfSubjects1}: The number of subjects by stage for
//'       the treatment group.
//'
//'     - \code{numberOfMilestone1}: The number of subjects reaching
//'       milestone by stage for the active treatment group.
//'
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{numberOfMilestone2}: The number of subjects reaching
//'       milestone by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfMilestone1}: The expected number of subjects
//'       reaching milestone for the active treatment group.
//'
//'     - \code{expectedNumberOfEvents2}: The expected number of events for
//'       control group.
//'
//'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
//'       for the control group.
//'
//'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
//'       for the control group.
//'
//'     - \code{expectedNumberOfMilestone2}: The expected number of subjects
//'       reaching milestone for the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' rmpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              rmstDiffLower = -2, rmstDiffUpper = 2,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List rmpowerequiv(const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.05,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const double milestone = NA_REAL,
                  const double rmstDiffLower = NA_REAL,
                  const double rmstDiffUpper = NA_REAL,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& lambda1 = NA_REAL,
                  const NumericVector& lambda2 = NA_REAL,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  const double accrualDuration = NA_REAL,
                  const double followupTime = NA_REAL,
                  const bool fixedFollowup = 0,
                  const NumericVector& spendingTime = NA_REAL,
                  const double studyDuration = NA_REAL) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(rmstDiffLower)) {
    stop("rmstDiffLower must be provided");
  }

  if (std::isnan(rmstDiffUpper)) {
    stop("rmstDiffUpper must be provided");
  }

  if (rmstDiffLower >= rmstDiffUpper) {
    stop("rmstDiffLower must be less than rmstDiffUpper");
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (lambda1.size() != 1 && lambda1.size() != nintervals &&
      lambda1.size() != nsi) {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() != 1 && lambda2.size() != nintervals &&
      lambda2.size() != nsi) {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
    stop("followupTime must be provided");
  }

  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }

  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      LogicalVector efficacyStopping1(kMax, 1);
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector efficacyP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
  }

  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, li, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));

  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  // obtain the timing of interim analysis
  NumericVector time(kMax);
  NumericVector u0(1, studyDuration1);
  DataFrame rm = rmstat(u0, milestone, allocationRatioPlanned,
                        accrualTime, accrualIntensity,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        accrualDuration, followupTime, fixedFollowup);

  double rmst1 = sum(NumericVector(rm[12]));
  double rmst2 = sum(NumericVector(rm[13]));
  double rmstDiff = sum(NumericVector(rm[14]));
  NumericVector theta(kMax, rmstDiff);
  double maxInformation = sum(NumericVector(rm[18]));

  // information at milestone
  u0[0] = milestone + 1.0e-6;
  rm = rmstat(u0, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double information1 = sum(NumericVector(rm[18]));
  if (informationRates1[0] <= information1/maxInformation) {
    std::string str1 = "The first information rate must exceed that";
    std::string str2 = "at the milestone time\n";
    std::string str3 = "The information at the milestone time is";
    std::string str4 = "The maximum information at study end is";
    std::string errmsg = str1 + " " + str2 +
      str3 + " " + std::to_string(information1) + "\n" +
      str4 + " " + std::to_string(maxInformation) + "\n";
    stop(errmsg);
  }

  NumericVector I = maxInformation*informationRates1;


  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup);
              return sum(NumericVector(rm[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    information1 = I[i];
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  };
  time[kMax-1] = studyDuration1;

  rm = rmstat(time, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double phi = allocationRatioPlanned/(allocationRatioPlanned+1);
  NumericVector nsubjects = NumericVector(rm[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
  NumericVector nevents = NumericVector(rm[2]);
  NumericVector nevents1 = NumericVector(rm[3]);
  NumericVector nevents2 = NumericVector(rm[4]);
  NumericVector ndropouts = NumericVector(rm[5]);
  NumericVector ndropouts1 = NumericVector(rm[6]);
  NumericVector ndropouts2 = NumericVector(rm[7]);
  NumericVector nmilestone = NumericVector(rm[9]);
  NumericVector nmilestone1 = NumericVector(rm[10]);
  NumericVector nmilestone2 = NumericVector(rm[11]);

  // calculate cumulative rejection probability under H1
  NumericVector theta10 = rep(rmstDiffLower, kMax);
  NumericVector theta20 = rep(rmstDiffUpper, kMax);
  NumericVector b = criticalValues1;
  NumericVector l = b + theta10*sqrt(I);
  NumericVector u = -b + theta20*sqrt(I);

  List probs1 = exitprobcpp(pmax(l, li), li, theta, I);
  List probs2 = exitprobcpp(ui, pmin(u, ui), theta, I);

  NumericVector cpl = cumsum(NumericVector(probs1[0]));
  NumericVector cpu = cumsum(NumericVector(probs2[1]));

  // identify the last look with l[k] >= u[k] if it exists
  IntegerVector k = which(l >= u);
  NumericVector cp(kMax);
  if (k.size() == 0) {
    cp = cpl + cpu - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], theta[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cp[i] = cpl[i] + cpu[i] - ca[i];
      } else {
        cp[i] = cpl[i] + cpu[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H1
  NumericVector q(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      q[i] = cp[i];
    } else if (i<kMax-1) {
      q[i] = cp[i] - cp[i-1];
    } else {
      q[i] = 1 - cp[i-1];
    }
  }

  NumericVector rejectPerStage(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      rejectPerStage[i] = cp[i];
    } else {
      rejectPerStage[i] = cp[i] - cp[i-1];
    }
  }

  NumericVector efficacyRmstDiffLower = theta10 + b/sqrt(I);
  NumericVector efficacyRmstDiffUpper = theta20 - b/sqrt(I);

  // calculate cumulative rejection under H10
  List probs2H10 = exitprobcpp(ui, pmin(u, ui), theta10, I);

  NumericVector cplH10 = cumAlphaSpent;
  NumericVector cpuH10 = cumsum(NumericVector(probs2H10[1]));

  NumericVector cpH10(kMax);
  if (k.size() == 0) {
    cpH10 = cplH10 + cpuH10 - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], theta10[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH10[i] = cplH10[i] + cpuH10[i] - ca[i];
      } else {
        cpH10[i] = cplH10[i] + cpuH10[i] - 1;
      }
    }
  }

  // calculate cumulative rejection under H20
  List probs1H20 = exitprobcpp(pmax(l, li), li, theta20, I);

  NumericVector cplH20 = cumsum(NumericVector(probs1H20[0]));
  NumericVector cpuH20 = cumAlphaSpent;

  NumericVector cpH20(kMax);
  if (k.size() == 0) {
    cpH20 = cplH20 + cpuH20 - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    NumericVector uH20 = -b;
    List a = exitprobcpp(l[idx], u[idx], theta20[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH20[i] = cplH20[i] + cpuH20[i] - ca[i];
      } else {
        cpH20[i] = cplH20[i] + cpuH20[i] - 1;
      }
    }
  }

  double overallReject = cp[kMax-1];
  double expectedNumberOfEvents = sum(q*nevents);
  double expectedNumberOfSubjects = sum(q*nsubjects);
  double expectedNumberOfEvents1 = sum(q*nevents1);
  double expectedNumberOfDropouts1 = sum(q*ndropouts1);
  double expectedNumberOfSubjects1 = sum(q*nsubjects1);
  double expectedNumberOfMilestone1 = sum(q*nmilestone1);
  double expectedNumberOfEvents2 = sum(q*nevents2);
  double expectedNumberOfDropouts2 = sum(q*ndropouts2);
  double expectedNumberOfSubjects2 = sum(q*nsubjects2);
  double expectedNumberOfMilestone2 = sum(q*nmilestone2);
  double expectedStudyDuration = sum(q*time);
  double expectedInformation = sum(q*I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["rmstDiffLower"] = rmstDiffLower,
    _["rmstDiffUpper"] = rmstDiffUpper,
    _["rmst1"] = rmst1,
    _["rmst2"] = rmst2,
    _["rmstDiff"] = rmstDiff,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup);

  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["rejectPerStage"] = rejectPerStage,
    _["cumulativeRejection"] = cp,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["cumulativeAttainedAlphaH10"] = cpH10,
    _["cumulativeAttainedAlphaH20"] = cpH20,
    _["numberOfEvents"] = nevents,
    _["numberOfDropouts"] = ndropouts,
    _["numberOfSubjects"] = nsubjects,
    _["numberOfMilestone"] = nmilestone,
    _["analysisTime"] = time,
    _["efficacyRmstDiffLower"] = efficacyRmstDiffLower,
    _["efficacyRmstDiffUpper"] = efficacyRmstDiffUpper,
    _["efficacyP"] = efficacyP,
    _["information"] = I);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["stratumFraction"] = stratumFraction,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["spendingTime"] = spendingTime);

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["numberOfMilestone1"] = nmilestone1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["numberOfMilestone2"] = nmilestone2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedNumberOfMilestone1"] = expectedNumberOfMilestone1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2,
    _["expectedNumberOfMilestone2"] = expectedNumberOfMilestone2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "rmpowerequiv";

  return result;
}


//' @title Sample Size for Equivalence in Restricted Mean Survival Time
//' Difference
//' @description Obtains the sample size for equivalence in restricted
//' mean survival time difference.
//'
//' @param beta The type II error.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffLower The lower equivalence limit of restricted mean
//'   survival time difference.
//' @param rmstDiffUpper The upper equivalence limit of restricted mean
//'   survival time difference.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return An S3 class \code{rmpowerequiv} object
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpowerequiv}}
//'
//' @examples
//'
//' rmsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   milestone = 18,
//'                   rmstDiffLower = -2, rmstDiffUpper = 2,
//'                   allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'                   accrualIntensity = 26/9*seq(1, 9),
//'                   piecewiseSurvivalTime = c(0, 6),
//'                   stratumFraction = c(0.2, 0.8),
//'                   lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'                   gamma1 = -log(1-0.05)/12,
//'                   gamma2 = -log(1-0.05)/12, accrualDuration = NA,
//'                   followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List rmsamplesizeequiv(const double beta = 0.2,
                       const int kMax = 1,
                       const NumericVector& informationRates = NA_REAL,
                       const NumericVector& criticalValues = NA_REAL,
                       const double alpha = 0.05,
                       const std::string typeAlphaSpending = "sfOF",
                       const double parameterAlphaSpending = NA_REAL,
                       const NumericVector& userAlphaSpending = NA_REAL,
                       const double milestone = NA_REAL,
                       const double rmstDiffLower = NA_REAL,
                       const double rmstDiffUpper = NA_REAL,
                       const double allocationRatioPlanned = 1,
                       const NumericVector& accrualTime = 0,
                       const NumericVector& accrualIntensity = NA_REAL,
                       const NumericVector& piecewiseSurvivalTime = 0,
                       const NumericVector& stratumFraction = 1,
                       const NumericVector& lambda1 = NA_REAL,
                       const NumericVector& lambda2 = NA_REAL,
                       const NumericVector& gamma1 = 0,
                       const NumericVector& gamma2 = 0,
                       double accrualDuration = NA_REAL,
                       double followupTime = NA_REAL,
                       const bool fixedFollowup = 0,
                       const NumericVector& interval =
                         NumericVector::create(0.001, 240),
                         const NumericVector& spendingTime = NA_REAL,
                         const bool rounding = 1) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi);


  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(rmstDiffLower)) {
    stop("rmstDiffLower must be provided");
  }

  if (std::isnan(rmstDiffUpper)) {
    stop("rmstDiffUpper must be provided");
  }

  if (rmstDiffLower >= rmstDiffUpper) {
    stop("rmstDiffLower must be less than rmstDiffUpper");

  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }

  if (!std::isnan(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }

  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }

    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }

  if (fixedFollowup && std::isnan(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (interval.size() != 2) {
    stop("interval must have 2 elements");
  }

  if (interval[0] < 0) {
    stop("lower limit of interval must be positive");
  }

  if (interval[0] >= interval[1]) {
    stop("upper limit must be greater than lower limit for interval");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      LogicalVector efficacyStopping1(kMax, 1);
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }


  std::string unknown;
  // search for the solution according to the input
  if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualDuration";
  } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
    unknown = "followupTime";
  } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }


  NumericVector b = criticalValues1;
  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);

  NumericVector rmsts1(nstrata), rmsts2(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  for (int h=0; h<nstrata; h++) {
    IntegerVector l = h*nintervals + l1;
    NumericVector lam1 = lambda1x[l];
    NumericVector lam2 = lambda2x[l];
    rmsts1[h] = rmst(0, milestone, piecewiseSurvivalTime, lam1);
    rmsts2[h] = rmst(0, milestone, piecewiseSurvivalTime, lam2);
  }

  double rmst1 = sum(stratumFraction*rmsts1);
  double rmst2 = sum(stratumFraction*rmsts2);
  double rmstDiff = rmst1 - rmst2;
  double theta10 = rmstDiffLower, theta20 = rmstDiffUpper;
  NumericVector theta(kMax, rmstDiff);

  List design = getDesignEquiv(
    beta, NA_REAL, theta10, theta20, rmstDiff,
    kMax, informationRates1, criticalValues1,
    alpha, asf, asfpar, userAlphaSpending, spendingTime1);

  DataFrame overallResults = DataFrame(design["overallResults"]);
  double maxInformation = overallResults["information"];
  double studyDuration;

  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            unknown, maxInformation](double aval)-> double{
              NumericVector accrualIntensity1 = clone(accrualIntensity);
              double dur1=0, dur2=0;

              if (unknown == "accrualDuration") {
                dur1 = aval;
                dur2 = followupTime;
              } else if (unknown == "followupTime") {
                dur1 = accrualDuration;
                dur2 = aval;
              } else if (unknown == "accrualIntensity") {
                dur1 = accrualDuration;
                dur2 = followupTime;
                accrualIntensity1 = aval*accrualIntensity;
              }

              // obtain the maximum information at study end
              NumericVector u0(1, dur1 + dur2);
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup);
              return sum(NumericVector(rm[18])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
    accrualDuration = brent(f, lower, interval[1], 1.0e-6);
    studyDuration = accrualDuration + followupTime;
  } else if (unknown == "followupTime") {
    if (f(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[18])) - maxInformation;
                };

      accrualDuration = brent(g, accrualDuration, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // adjust follow-up time
      double lower =  std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(f, lower, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  } else if (unknown == "accrualIntensity") {
    if (accrualDuration + followupTime <= milestone) {
      std::string str1 = "NOTE: The followupTime is too short for ";
      std::string str2 = "anyone to reach milestone.";
      std::string str3 = "NOTE: followupTime is increased to ";
      std::string str4 = "milestone.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";
      followupTime = milestone;
    }
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
    studyDuration = accrualDuration + followupTime;
  }


  // output the results
  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup);
                    return sum(NumericVector(rm[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // the follow-up time cannot be adjusted for fixed follow-up
        // adjust study duration to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup);
                    return sum(NumericVector(rm[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  List result = rmpowerequiv(
    kMax, informationRates1, criticalValues1,
    alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, milestone, rmstDiffLower, rmstDiffUpper,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  return result;
}
