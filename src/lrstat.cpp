#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Number of enrolled subjects
//' @description Obtains the number of subjects enrolled by given calendar
//' times.
//'
//' @param time A vector of calendar times at which to calculate the number
//' of enrolled subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_accrualDuration
//'
//' @return A vector of total number of subjects enrolled by the
//' specified calendar times.
//'
//' @examples
//' # Example 1: Uniform enrollment with 20 patients per month for 12 months.
//' accrual(time = 3, accrualTime = 0, accrualIntensity = 20,
//'         accrualDuration = 12)
//'
//' # Example 2: Piecewise accrual, 10 patients per month for the first
//' # 3 months, and 20 patients per month thereafter. Patient recruitment
//' # ends at 12 months for the study.
//' accrual(time = c(2, 9), accrualTime = c(0, 3),
//'         accrualIntensity = c(10, 20), accrualDuration = 12)
//'
//' @export
// [[Rcpp::export]]
NumericVector accrual(const NumericVector& time = NA_REAL,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const double accrualDuration = NA_REAL) {

  int i, j, k = time.size();
  NumericVector n(k);

  // up to end of enrollment
  NumericVector t = pmax(pmin(time, accrualDuration), 0.0);

  // identify the time interval containing t
  IntegerVector m = pmax(findInterval2(t, accrualTime), 1);

  // sum up patients enrolled in each interval up to t
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        n[i] += accrualIntensity[j]*(accrualTime[j+1] - accrualTime[j]);
      } else {
        n[i] += accrualIntensity[j]*(t[i] - accrualTime[j]);
      }
    }
  }

  return n;
}

//' @title Probability of being at risk
//' @description Obtains the probability of being at risk at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//' probability of being at risk.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of being at risk at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//' patrisk(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'         lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector patrisk(const NumericVector& time = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda = NA_REAL,
                      const NumericVector& gamma = 0) {

  // identify the time interval containing the specified analysis time
  IntegerVector m = pmax(findInterval2(time, piecewiseSurvivalTime), 1);

  int i, j, k = time.size(), J = lambda.size();

  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }

  NumericVector t = piecewiseSurvivalTime;

  // sum up cumulative hazard up to time
  NumericVector a(k);
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        a[i] += lg[j]*(t[j+1] - t[j]);
      } else {
        a[i] += lg[j]*(time[i] - t[j]);
      }
    }
  }

  return exp(-a);
}


//' @title Probability of having an event
//' @description Obtains the probability of having an event at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//' probability of having an event.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of having an event at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//' pevent(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'        lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector pevent(const NumericVector& time = NA_REAL,
                     const NumericVector& piecewiseSurvivalTime = 0,
                     const NumericVector& lambda = NA_REAL,
                     const NumericVector& gamma = 0) {

  // identify the time interval containing the specified analysis time
  IntegerVector m = pmax(findInterval2(time, piecewiseSurvivalTime), 1);

  int i, j, k = time.size(), J = lambda.size();

  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }

  // sum up cumulative hazard up to time
  NumericVector t = piecewiseSurvivalTime;
  NumericVector n = patrisk(t, t, lambda, gamma);
  NumericVector a(k);
  double p;

  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        p = lambda[j]/lg[j]*(1 - exp(-lg[j]*(t[j+1] - t[j])));
      } else {
        p = lambda[j]/lg[j]*(1 - exp(-lg[j]*(time[i] - t[j])));
      }
      a[i] += n[j]*p;
    }
  }

  return a;
}


//' @title Integrated event probability over an interval with constant hazard
//' @description Obtains the integration probability of having an event
//' during an interval with constant hazard.
//'
//' @param j The analysis time interval with constant hazard.
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval with constant hazard for a treatment
//' group with specified piecewise exponential survival and dropout
//' distributions.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//' hd(j = 1, t1 = 1, t2 = 3, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double hd(const int j = NA_INTEGER,
          const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {

  int j1 = j-1;

  // lower bound of time interval j for the piecewise exponential distribution
  NumericVector t0 = NumericVector::create(piecewiseSurvivalTime[j1]);

  // probability of being at risk at the start of interval j
  NumericVector n0 = patrisk(t0, piecewiseSurvivalTime, lambda, gamma);

  // probability of having an event at the start of interval j
  NumericVector d0 = pevent(t0, piecewiseSurvivalTime, lambda, gamma);


  int J = lambda.size();

  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }

  // integration of conditional probability of having an event over (t1,t2)
  // given survival at the start of interval j
  double q1 = (exp(-lg[j1]*(t1-t0[0])) - exp(-lg[j1]*(t2-t0[0])))/lg[j1];
  double q = lambda[j1]/lg[j1] * (t2-t1 - q1);

  // sum up the integration for the already failed and to-be-failed
  return d0[0]*(t2-t1) + n0[0]*q;
}


//' @title Integrated event probability over an interval
//' @description Obtains the integration of the probability of having an event
//' during an interval. The specified analysis time interval can span more
//' than one analysis time interval with constant hazard.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval for a treatment group with specified
//' piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//' pd(t1 = 1, t2 = 8, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double pd(const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {

  // identify the analysis time intervals containing t1 and t2
  NumericVector t12 = NumericVector::create(t1, t2);
  IntegerVector j12 = pmax(findInterval2(t12, piecewiseSurvivalTime), 1) - 1;

  NumericVector t = piecewiseSurvivalTime;

  int j, j1=j12[0], j2=j12[1];

  // sum up the integrated event probabilities across analysis time intervals
  double a=0, x;
  for (j=j1; j<=j2; j++) {
    if (j1==j2) {
      x = hd(j+1, t1, t2, t, lambda, gamma);
    } else if (j==j1) {
      x = hd(j+1, t1, t[j+1], t, lambda, gamma);
    } else if (j==j2) {
      x = hd(j+1, t[j], t2, t, lambda, gamma);
    } else {
      x = hd(j+1, t[j], t[j+1], t, lambda, gamma);
    }
    a += x;
  }

  return a;
}


//' @title Number of patients enrolled during an interval and having an event
//' by specified calendar times
//' @description Obtains the number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//' of patients having an event.
//' @param u1 Lower bound of the accrual time interval.
//' @param u2 Upper bound of the accrual time interval.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times for a given treatment group had the enrollment being
//' restricted to the treatment group. By definition, we must have
//' \code{time >= u2}.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise accrual, 10 patients per month for the first 3 months, and
//' # 20 patients per month thereafter. Piecewise exponential survival with
//' # hazard 0.0533 in the first 6 months, and hazard 0.0309 thereafter,
//' # and 5% dropout by the end of 1 year.
//' ad(time = c(9, 15), u1 = 1, u2 = 8, accrualTime = c(0, 3),
//'    accrualIntensity = c(10, 20), piecewiseSurvivalTime=c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector ad(const NumericVector& time = NA_REAL,
                 const double u1 = NA_REAL,
                 const double u2 = NA_REAL,
                 const NumericVector& accrualTime = 0,
                 const NumericVector& accrualIntensity = NA_REAL,
                 const NumericVector& piecewiseSurvivalTime = 0,
                 const NumericVector& lambda = NA_REAL,
                 const NumericVector& gamma = 0) {

  // identify the accrual time intervals containing u1 and u2
  NumericVector u12 = NumericVector::create(u1, u2);
  IntegerVector j12 = pmax(findInterval2(u12, accrualTime), 1) - 1;

  NumericVector u = accrualTime;

  int i, j, j1=j12[0], j2=j12[1], k=time.size();

  NumericVector a(k);

  // sum up the number of patients with event across accrual time intervals
  double t, x;
  for (i=0; i<k; i++) {
    t = time[i];
    for (j=j1; j<=j2; j++) {
      if (j1==j2) {
        x = pd(t-u2, t-u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j==j1) {
        x = pd(t-u[j+1], t-u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j==j2) {
        x = pd(t-u2, t-u[j], piecewiseSurvivalTime, lambda, gamma);
      } else {
        x = pd(t-u[j+1], t-u[j], piecewiseSurvivalTime, lambda, gamma);
      }
      a[i] += accrualIntensity[j]*x;
    }
  }

  return a;
}


//' @title Number of subjects at risk
//' @description Obtains the number of subjects at risk at given analysis
//' times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//' of patients at risk.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients at risk at the specified
//' analysis times for each treatment group.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' natrisk(time = c(9, 24), allocationRatioPlanned = 1, accrualTime = c(0, 3),
//'         accrualIntensity = c(10, 20), piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix natrisk(const NumericVector& time = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double minFollowupTime = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {

  // truncate the analysis time by the maximum follow-up
  NumericVector t = pmin(time, maxFollowupTime);

  // enrollment time
  NumericVector u = pmin(accrualDuration+minFollowupTime-t, accrualDuration);

  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  // number of patients at risk in each treatment group
  int k = time.size();
  NumericMatrix n(k, 2);
  n(_, 0) = phi*a*patrisk(t, piecewiseSurvivalTime, lambda1, gamma1);
  n(_, 1) = (1-phi)*a*patrisk(t, piecewiseSurvivalTime, lambda2, gamma2);

  return n;
}


//' @title Number of subjects having an event
//' @description Obtains the number of subjects having an event by given
//' analysis times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//' of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the specified
//' analysis times for each treatment group.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' nevent(time = c(9, 24), allocationRatioPlanned = 1, accrualTime = c(0, 3),
//'        accrualIntensity = c(10, 20), piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 12, minFollowupTime = 18, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix nevent(const NumericVector& time = NA_REAL,
                     const double allocationRatioPlanned = 1,
                     const NumericVector& accrualTime = 0,
                     const NumericVector& accrualIntensity = NA_REAL,
                     const NumericVector& piecewiseSurvivalTime = 0,
                     const NumericVector& lambda1 = NA_REAL,
                     const NumericVector& lambda2 = NA_REAL,
                     const NumericVector& gamma1 = 0,
                     const NumericVector& gamma2 = 0,
                     const double accrualDuration = NA_REAL,
                     const double minFollowupTime = NA_REAL,
                     const double maxFollowupTime = NA_REAL) {

  // truncate the analysis time by the maximum follow-up
  NumericVector t = pmin(time, maxFollowupTime);

  // enrollment time
  NumericVector u = pmin(accrualDuration+minFollowupTime-t, accrualDuration);

  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  // number of patients having an event in each treatment group
  NumericVector u1(1);
  u1[0] = accrualDuration + minFollowupTime;

  int i, k = time.size();
  NumericMatrix d(k, 2);

  NumericVector d1(k), d2(k);
  d1 = a*pevent(t, piecewiseSurvivalTime, lambda1, gamma1);
  d2 = a*pevent(t, piecewiseSurvivalTime, lambda2, gamma2);

  for (i=0; i<k; i++) {
    d(i,0) = phi*(d1[i] + ad(u1, u[i], accrualDuration, accrualTime,
                  accrualIntensity, piecewiseSurvivalTime, lambda1, gamma1)[0]);
    d(i,1) = (1-phi)*(d2[i] + ad(u1, u[i], accrualDuration, accrualTime,
              accrualIntensity, piecewiseSurvivalTime, lambda2, gamma2)[0]);
  }

  return d;
}


//' @title Number of subjects having an event by calendar time
//' @description Obtains the number of subjects having an event by given
//' calendar times for each treatment group.
//'
//' @param time A vector of calendar times at which to calculate the number
//' of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the specified
//' calendar times for each treatment group.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' nevent2(time = c(9, 24), allocationRatioPlanned = 1, accrualTime = c(0, 3),
//'         accrualIntensity = c(10, 20), piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix nevent2(const NumericVector& time = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double minFollowupTime = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {

  // truncate the calendar time by study end
  NumericVector t = pmin(time, accrualDuration + minFollowupTime);

  // enrollment time
  NumericVector u = pmin(pmax(t - maxFollowupTime, 0.0), accrualDuration);
  NumericVector w = pmin(t, accrualDuration);

  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  // number of patients having an event in each treatment group
  NumericVector s(1), v(1);
  s[0] = maxFollowupTime;

  int i, k = time.size();
  NumericMatrix d(k, 2);

  NumericVector d1(k), d2(k);
  d1 = a*pevent(s, piecewiseSurvivalTime, lambda1, gamma1)[0];
  d2 = a*pevent(s, piecewiseSurvivalTime, lambda2, gamma2)[0];

  for (i=0; i<k; i++) {
    v[0] = t[i];
    d(i,0) = phi*(d1[i] + ad(v, u[i], w[i], accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda1, gamma1)[0]);
    d(i,1) = (1-phi)*(d2[i] + ad(v, u[i], w[i], accrualTime, accrualIntensity,
              piecewiseSurvivalTime, lambda2, gamma2)[0]);
  }

  return d;
}

//' @title Number of subjects having an event and log-rank statistics
//' @description Obtains the number of subjects having an event in each
//' treatment group, the mean and variance of weighted log-rank test score
//' statistic at given calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of events and the mean and variance of log-rank test score statistic.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//' @param predictEventOnly Whether to predict the number of events only.
//'   Defaults to 0 for obtaining log-rank test score statistic mean
//'   and variance.
//'
//' @return A data frame of the number of patients enrolled, the number of
//' patients having an event overall and in each treatment group, the mean and
//' variance of weighted log-rank test score statistic at the specified
//' calendar times by stratum.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' lrstat(time = c(22, 40), allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 9),
//'        accrualIntensity = c(26/9*seq(1, 9), 26),
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
DataFrame lrstat(const NumericVector& time = NA_REAL,
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
                 const double rho1 = 0,
                 const double rho2 = 0,
                 const int numSubintervals = 300,
                 const bool predictEventOnly = 0) {

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();

  if (nintervals*nstrata != lambda1.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda1");
  }

  if (nintervals*nstrata != lambda2.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("gamma1 must be a scalar or have the same length as lambda1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
    stop("gamma2 must be a scalar or have the same length as lambda2");
  }


  double minFollowupTime = followupTime;
  double maxFollowupTime;

  // obtain the follow-up time for the first enrolled subject
  if (fixedFollowup) {
    maxFollowupTime = minFollowupTime;
  } else {
    maxFollowupTime = accrualDuration + minFollowupTime;
  }

  int i, j, j1, k = time.size(), h;
  int sk = nstrata*k;

  IntegerVector stratum(sk);
  NumericVector times(sk);
  NumericVector nsubjects(sk);
  NumericMatrix nevents(sk, 2);
  NumericVector nevents1(sk), nevents2(sk), neventst(sk);
  NumericVector uscore(sk), vscore(sk);

  double frac, accrualDuration0, minFollowupTime0, maxFollowupTime0, inc;
  IntegerVector l(nintervals);
  IntegerVector l1 = Range(0, nintervals-1);
  NumericVector lam1(nintervals), lam2(nintervals);
  IntegerVector q = Range(0, numSubintervals);
  NumericVector q1 = as<NumericVector>(q);
  Range q2 = Range(0, numSubintervals-1), c0 = Range(0,0), c1 = Range(1,1);
  NumericVector t(numSubintervals+1);
  NumericMatrix xatrisk(numSubintervals+1, 2);
  NumericMatrix xevent(numSubintervals+1, 2);
  NumericVector atrisk1(numSubintervals), atrisk2(numSubintervals),
  atriskt(numSubintervals), event1(numSubintervals), event2(numSubintervals),
  eventt(numSubintervals), km(numSubintervals), w(numSubintervals);

  NumericVector s = pmin(time, accrualDuration + minFollowupTime);
  NumericVector a = accrual(s, accrualTime, accrualIntensity, accrualDuration);

  NumericMatrix x(k, 2);

  DataFrame df;


  for (h=0; h<nstrata; h++) {
    frac = stratumFraction[h];
    l = h*nintervals + l1;
    lam1 = lambda1[l];
    lam2 = lambda2[l];

    // number of events in the stratum at the specified calendar times
    x = nevent2(s, allocationRatioPlanned, accrualTime, frac*accrualIntensity,
                piecewiseSurvivalTime, lam1, lam2, gamma1, gamma2,
                accrualDuration, minFollowupTime, maxFollowupTime);

    // loop over each calendar time point of interest
    for (j=0; j<k; j++) {
      j1 = h*k + j;

      // obtain number of enrolled subjects and subjects having an event
      nsubjects[j1] = frac*a[j];
      nevents(j1, _) = x.row(j);

      // approximate the mean and variance of weighted log-rank test
      // score statistic
      if (!predictEventOnly) {

        // modify the study design at the calendar time of interest
        accrualDuration0 = std::min(s[j], accrualDuration);
        minFollowupTime0 = std::max(s[j] - accrualDuration, 0.0);
        maxFollowupTime0 = std::min(s[j], maxFollowupTime);

        // partition the follow-up period into small sub-intervals
        inc = maxFollowupTime0/numSubintervals;
        t = q1*inc;

        // obtain number of patients at risk and the number of patients having
        // an event at each analysis time point
        xatrisk = natrisk(t, allocationRatioPlanned,
                          accrualTime, frac*accrualIntensity,
                          piecewiseSurvivalTime, lam1, lam2, gamma1, gamma2,
                          accrualDuration0, minFollowupTime0, maxFollowupTime0);

        xevent = nevent(t, allocationRatioPlanned,
                        accrualTime, frac*accrualIntensity,
                        piecewiseSurvivalTime, lam1, lam2, gamma1, gamma2,
                        accrualDuration0, minFollowupTime0, maxFollowupTime0);

        // number of patients at risk at start of each analysis time interval
        atrisk1 = xatrisk(q2, c0);
        atrisk2 = xatrisk(q2, c1);
        atriskt = atrisk1 + atrisk2;

        // number of patients having an event in each analysis time interval
        event1 = diff(xevent(_, 0));
        event2 = diff(xevent(_, 1));
        eventt = event1 + event2;

        // Kaplan-Meier estimates of survival probabilities at the start of
        // each analysis time interval
        km[0] = 1;
        for (i=1; i<numSubintervals; i++) {
          km[i] = km[i-1]*(1 - eventt[i-1]/atriskt[i-1]);
        }

        // vector of Fleming-Harrington weights
        w = pow(km,rho1)*pow(1-km,rho2);

        // mean of the weighted log-rank test score statistic
        uscore[j1] = sum(w * (event1 - eventt*atrisk1/atriskt));

        // variance of the weighted log-rank test score statistic
        vscore[j1] = sum(w*w * eventt*(atrisk1*atrisk2/pow(atriskt,2)));
      }
    }
  }

  // number of subjects having an event in each treatment group and overall
  nevents1 = nevents(_, 0);
  nevents2 = nevents(_, 1);
  neventst = nevents1 + nevents2;

  // stratum and time
  for (h=0; h<nstrata; h++) {
    for (j=0; j<k; j++) {
      stratum[h*k+j] = h+1;
      times[h*k+j] = time[j];
    }
  }


  // output the requested information
  if (predictEventOnly) {
    df = DataFrame::create(_["stratum"] = stratum,
                           _["time"] = times,
                           _["subjects"] = nsubjects,
                           _["nevents"] = neventst,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2);
  } else {
    df = DataFrame::create(_["stratum"] = stratum,
                           _["time"] = times,
                           _["subjects"] = nsubjects,
                           _["nevents"] = neventst,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore);
  }

  return df;
}


//' @title Calendar times for target number of events
//' @description Obtains the calendar times to reach the target number of
//' subjects having an event.
//'
//' @param nevents A vector of target number of events.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//'
//' @return A vector of calendar times expected to yield the target
//' number of events.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' caltime(nevents = c(24, 80), allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 9),
//'         accrualIntensity = c(26/9*seq(1, 9), 26),
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
NumericVector caltime(const NumericVector& nevents = NA_REAL,
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

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();

  if (nintervals*nstrata != lambda1.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda1");
  }

  if (nintervals*nstrata != lambda2.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("gamma1 must be a scalar or have the same length as lambda1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
    stop("gamma2 must be a scalar or have the same length as lambda2");
  }

  NumericVector t0(1);  // convert a number into a vector
  DataFrame lr;
  double event;

  // Lambda function
  auto f = [&](double t) {
    t0[0] = t;
    lr = lrstat(t0, allocationRatioPlanned, accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                0, 0, 300, 1);
    return sum(NumericVector(lr[3])) - event;  // sum events across strata;
  };

  int i, k = nevents.size();
  double studyTime = accrualDuration + followupTime;
  NumericVector time(k);

  event = max(nevents);
  if (f(studyTime) < 0.0) {
    stop("followupTime is too short to reach the target number of events");
  }

  for (i=0; i<k; i++) {
    // match the predicted number of events to the target
    event = std::max(nevents[i], 0.0);
    time[i] = brent(f, 0.0001, studyTime, 0.0001);
  }

  return time;
}


#include <Rcpp.h>
using namespace Rcpp;

//' @title Stagewise exit probabilities
//' @description Obtains the stagewise exit probabilities for both efficacy and
//' futility stopping.
//'
//' @param b Upper boundaries on the z-test statistic scale.
//' @param a Lower boundaries on the z-test statistic scale. Defaults to
//' \code{c(rep(-Inf, kMax-1), b[kMax])} if left unspecified, where
//' \code{kMax = length(b)}.
//' @param theta Stagewise parameter of interest, e.g., \code{-U/V} for
//' weighted log-rank test, where \code{U} is the mean and \code{V} is
//' the variance of the weighted log-rank test score statistic at each stage.
//' For proportional hazards and conventional log-rank test, use the
//' scalar input, \code{theta = -log(HR)}.
//' @param I Stagewise cumulative information, e.g., \code{V}, the variance
//' of the weighted log-rank test score statistic at each stage. For
//' conventional log-rank test, information can be approximated by
//' \code{phi*(1-phi)*D}, where \code{phi} is the probability of being
//' allocated to the active arm, and \code{D} is the total number of events at
//' each stage.
//' @param r Integer value controlling grid for numerical integration as in
//' Jennison and Turnbull (2000). Defaults to 18. Specify a larger number
//' for greater accuracy.
//'
//' @return A list of stagewise exit probabilities: one vector for efficacy
//' stopping probabilities, and the other vector for futility stopping
//' probabilities.
//'
//' @examples
//' exitprob(b = c(3.471, 2.454, 2.004), theta = -log(0.6),
//'          I = c(50, 100, 150)/4)
//'
//' @export
// [[Rcpp::export]]
List exitprob(const NumericVector& b = NA_REAL,
              NumericVector a = NA_REAL,
              NumericVector theta = NA_REAL,
              const NumericVector& I = NA_REAL,
              const int r = 18) {

  // variable declarations
  // kMax is the total number of stages
  // m0, z0, h0 for the previous stage
  // m, z, h for the current stage
  int kMax=b.size(), r1=6*r-1, r2=12*r-3, i0, i1=0, i2=r1-1, i, j,
    m0=r2, m1=r1, m=r2;
  double t, tlower, tupper, xlower, xupper;

  NumericVector sqrtI(kMax), thetaSqrtI(kMax), thetaI(kMax), dI(kMax),
  dThetaI(kMax), exitProbUpper(kMax), exitProbLower(kMax),
  shift(r1), x1(r1), x(r1), z0(r2), z(r2), w(r2), h0(r2), h(r2);

  // set default parameter values
  if (any(is_na(a))) {
    NumericVector tem(kMax);
    for (i=0; i<kMax; i++) {
      if (i<kMax-1) {
        tem[i] = R_NegInf;
      } else {
        tem[i] = b[i];
      }
    }
    a = tem;
  }


  if (theta.size()==1) {
    theta = rep(theta, kMax);
  }

  if (a.size()!=kMax || theta.size()!=kMax || I.size()!=kMax) {
    stop("The input parameters must have the same length");
  }


  // constant shifts relative to the means, use floating point computation
  for (i=0; i<r1; i++) {
    if (i < r-1) {
      shift[i] = -3 - 4*log(r/(i+1.0));
    } else if (i < 5*r) {
      shift[i] = -3 + 3*(i+1.0-r)/(2*r);
    } else {
      shift[i] = 3 + 4*log(r/(6*r-i-1.0));
    }
  }


  // obtain various vectors associated with theta and I
  for (j=0; j<kMax; j++) {
    sqrtI[j] = sqrt(I[j]);
    thetaSqrtI[j] = theta[j]*sqrtI[j];
    thetaI[j] = theta[j]*I[j];
    if (j==0) {
      dI[j] = I[j];
      dThetaI[j] = thetaI[j];
    } else {
      dI[j] = I[j] - I[j-1];
      dThetaI[j] = thetaI[j] - thetaI[j-1];
    }
  }

  // loop over stages
  for (j=0; j<kMax; j++) {

    // initialize x values
    for (i=0; i<r1; i++) {
      x1[i] = thetaSqrtI[j] + shift[i];
    }

    // trim off x values outside (a[j], b[j])
    // trim from below
    if (a[j] >= x1[0]) {
      i1 = 0;
      while (x1[i1] <= a[j]) {
        i1++;
      }
      i1--;
      xlower = a[j];  // lower bound on x
    } else {
      i1 = 0;
      xlower = x1[0];
    }

    // trim from above
    if (b[j] <= x1[r1-1]) {
      i2 = r1-1;
      while (x1[i2] >= b[j]) {
        i2--;
      }
      i2++;
      xupper = b[j];  // upper bound on x
    } else {
      i2 = r1-1;
      xupper = x1[r1-1];
    }

    // save the trimmed portion to x
    m1 = i2 - i1 + 1;
    x[0] = xlower;
    x[m1-1] = xupper;
    for (i=1; i<m1-1; i++) {
      x[i] = x1[i+i1];
    }

    // derive the grid points for z
    m = 2*m1 - 1;

    // odd grid points;
    for (i=0; i<m1; i++) {
      z[2*i] = x[i];
    }

    // even grid points;
    for (i=0; i<m1-1; i++) {
      z[2*i+1] = (z[2*i] + z[2*i+2])/2;
    }


    // derive the weights
    w[0] = 1.0/6*(z[2] - z[0]);

    for (i0=1; i0<=m1-2; i0++) {
      i = 2*i0;
      w[i] = 1.0/6*(z[i+2] - z[i-2]);
    }

    for (i0=1; i0<=m1-1; i0++) {
      i = 2*i0-1;
      w[i] = 4.0/6*(z[i+1] - z[i-1]);
    }

    w[m-1] = 1.0/6*(z[m-1] - z[m-3]);


    // first stage is easy
    if (j==0) {
      // exit probabilities
      exitProbUpper[j] = R::pnorm(-b[j] + thetaSqrtI[j], 0.0, 1.0, 1, 0);
      exitProbLower[j] = R::pnorm(a[j] - thetaSqrtI[j], 0.0, 1.0, 1, 0);

      // prepare h0, m0, z0 for the next stage
      if (kMax > 1) {
        for (i=0; i<m; i++) {
          h0[i] = w[i]*R::dnorm(z[i] - thetaSqrtI[j], 0.0, 1.0, 0);
        }

        m0 = m;
        z0 = z+0.0; // adding 0.0 to avoid passing by reference
      }

    } else {
      // calculate exit probabilities using h0 from the previous stage
      for (i0=0; i0<m0; i0++) {
        tupper = (z0[i0]*sqrtI[j-1] - b[j]*sqrtI[j] + dThetaI[j])/sqrt(dI[j]);
        tlower = (-z0[i0]*sqrtI[j-1] + a[j]*sqrtI[j] - dThetaI[j])/sqrt(dI[j]);
        exitProbUpper[j] += h0[i0]*R::pnorm(tupper, 0.0, 1.0, 1, 0);
        exitProbLower[j] += h0[i0]*R::pnorm(tlower, 0.0, 1.0, 1, 0);
      }

      // prepare h0, m0, z0 for the next stage
      if (j < kMax-1) {
        for (i=0; i<m; i++) {
          h[i] = 0;
          for (i0=0; i0<m0; i0++) {
            t = (z[i]*sqrtI[j] - z0[i0]*sqrtI[j-1] - dThetaI[j])/sqrt(dI[j]);
            h[i] += h0[i0]*R::dnorm(t, 0.0, 1.0, 0);
          }
          h[i] *= w[i]*sqrt(I[j]/dI[j]);  // factors invariant to i0
        }

        h0 = h+0.0; // adding 0.0 to avoid passing by reference
        m0 = m;
        z0 = z+0.0;
      }
    }

  }

  // return a list of stagewise exit probabilities
  return List::create(Named("exitProbUpper") = exitProbUpper,
                      Named("exitProbLower") = exitProbLower);

}


//' @title Log-rank test power
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a two-group survival design.
//'
//' @inheritParams param_kMax
//' @inheritParams param_informationRates
//' @inheritParams param_criticalValues
//' @inheritParams param_futilityBounds
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//'
//' @return A list of the overall and stagewise rejection probabilities, the
//' futility stoppig probabilities, the overall and stagewise expected number
//' of events, number of patients, and analysis time, the input accrual and
//' follow-up durations, and whether a fixed follow-up is used.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrpower(kMax = 2, informationRates = c(0.8, 1),
//'         criticalValues = c(2.250, 2.025),
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 9),
//'         accrualIntensity = c(26/9*seq(1, 9), 26),
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
List lrpower(const int kMax = NA_INTEGER,
             NumericVector informationRates = NA_REAL,
             const NumericVector& criticalValues = NA_REAL,
             NumericVector futilityBounds = NA_REAL,
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
             const double rho1 = 0,
             const double rho2 = 0,
             const int numSubintervals = 300) {

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();

  if (nintervals*nstrata != lambda1.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda1");
  }

  if (nintervals*nstrata != lambda2.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("gamma1 must be a scalar or have the same length as lambda1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
    stop("gamma2 must be a scalar or have the same length as lambda2");
  }


  // set default parameter values
  if (any(is_na(informationRates))) {
    IntegerVector tem = seq_len(kMax);
    informationRates = as<NumericVector>(tem)/(kMax+0.0);
  }

  if (kMax > 1) {
    if (any(is_na(futilityBounds))) {
      futilityBounds = rep(R_NegInf, kMax-1);
    }
  }

  int i, j, k;

  NumericVector u = criticalValues;
  NumericVector l(kMax);
  if (kMax > 1) {
    for (i=0; i<kMax; i++) {
      if (i<kMax-1) {
        l[i] = futilityBounds[i];
      } else {
        l[i] = u[i];
      }
    }
  } else {
    l = u;
  }

  NumericVector u0(1);
  DataFrame lr;
  NumericVector e0(kMax), time(kMax);

  // obtain the total number of events at study end
  u0[0] = accrualDuration + followupTime;
  lr = lrstat(u0, allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 1);


  // obtain the timing of interim analysis
  e0 = sum(NumericVector(lr[3]))*informationRates;
  time = caltime(e0, allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 accrualDuration, followupTime, fixedFollowup);

  // obtain the mean and variance of log-rank test score statistic at
  // each stage
  lr = lrstat(time, allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 0);


  // add up the mean and variance across strata;
  NumericVector nsubjects(kMax), nevents(kMax), uscore(kMax), vscore(kMax);
  for (i=0; i<nstrata; i++) {
    for (j=0; j<kMax; j++) {
      k = i*kMax+j;
      nsubjects[j] += NumericVector(lr[2])[k];
      nevents[j] += NumericVector(lr[3])[k];
      uscore[j] += NumericVector(lr[6])[k];
      vscore[j] += NumericVector(lr[7])[k];
    }
  }

  // compute the stagewise exit probabilities for efficacy and futility
  NumericVector theta = -uscore/vscore;
  List probs = exitprob(u, l, theta, vscore);

  // stagewise total exit probabilities
  NumericVector pu(kMax), pl(kMax), ptotal(kMax);
  pu = NumericVector(probs[0]);
  pl = NumericVector(probs[1]);
  ptotal = pu + pl;

  double overallReject = sum(pu);
  double expectedNumberOfEvents = sum(ptotal*nevents);
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedStudyDuration = sum(ptotal*time);

  List z;

  if (kMax > 1) {
    z = List::create(_["overallReject"] = overallReject,
                     _["rejectPerStage"] = pu,
                     _["futilityPerStage"] = pl,
                     _["eventsPerStage"] = nevents,
                     _["numberOfSubjects"] = nsubjects,
                     _["analysisTime"] = time,
                     _["expectedNumberOfEvents"] = expectedNumberOfEvents,
                     _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
                     _["expectedStudyDuration"] = expectedStudyDuration,
                     _["accrualDuration"] = accrualDuration,
                     _["followupTime"] = followupTime,
                     _["fixedFollowup"] = fixedFollowup);
  } else {
    z = List::create(_["overallReject"] = overallReject,
                     _["numberOfEvents"] = nevents,
                     _["numberOfSubjects"] = nsubjects,
                     _["analysisTime"] = time,
                     _["accrualDuration"] = accrualDuration,
                     _["followupTime"] = followupTime,
                     _["fixedFollowup"] = fixedFollowup);
  }

  return z;
}


//' @title Log-rank test sample size
//' @description Obtains the needed accrual duration given power and
//' follow-up time, or the needed follow-up time given power and
//' accrual duration in a two-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @inheritParams param_informationRates
//' @inheritParams param_criticalValues
//' @inheritParams param_futilityBounds
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//' @param interval The interval to search for the solution of
//' accrualDuration or followupDuration. Defaults to \code{c(0.001, 240)}.
//' Adjustment may be needed for non-monotone relationship with study power.
//'
//' @return A list of the overall and stagewise rejection probabilities,  the
//' futility stoppig probabilities, the overall and stagewise expected number
//' of events, number of patients, and analysis time, the input or calculated
//' accrual and follow-up durations, and whether a fixed follow-up is used.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' # Example 1: Obtains accrual duration given power and follow-up duration
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              criticalValues = c(2.250, 2.025),
//'              accrualTime = seq(0, 9),
//'              accrualIntensity = c(26/9*seq(1, 9), 26),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = NA,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains follow-up duration given power and accrual duration
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              criticalValues = c(2.250, 2.025),
//'              accrualTime = seq(0, 9),
//'              accrualIntensity = c(26/9*seq(1, 9), 26),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = 22,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List lrsamplesize(const double beta = 0.2,
                  const int kMax = NA_INTEGER,
                  NumericVector informationRates = NA_REAL,
                  const NumericVector& criticalValues = NA_REAL,
                  NumericVector futilityBounds = NA_REAL,
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
                  const double rho1 = 0,
                  const double rho2 = 0,
                  const int numSubintervals = 300,
                  const NumericVector& interval =
                    NumericVector::create(0.001, 240)) {

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();

  if (nintervals*nstrata != lambda1.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda1");
  }

  if (nintervals*nstrata != lambda2.size()) {
    stop("piecewiseSurvivalTime must have the same length as lambda2");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("gamma1 must be a scalar or have the same length as lambda1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
    stop("gamma2 must be a scalar or have the same length as lambda2");
  }


  // set default parameter values
  if (any(is_na(informationRates))) {
    IntegerVector tem = seq_len(kMax);
    informationRates = as<NumericVector>(tem)/kMax;
  }

  if (kMax > 1) {
    if (any(is_na(futilityBounds))) {
      futilityBounds = rep(R_NegInf, kMax-1);
    }
  }


  List lr;

  // search for the solution according to the input
  if (R_isnancpp(accrualDuration) && !R_isnancpp(followupTime)) {
    // Lambda function
    auto f = [&](double t) {
      lr = lrpower(kMax, informationRates, criticalValues, futilityBounds,
                   allocationRatioPlanned, accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, stratumFraction,
                   lambda1, lambda2, gamma1, gamma2,
                   t, followupTime, fixedFollowup,
                   rho1, rho2, numSubintervals);
      return double(lr[0]) - (1-beta);
    };
    accrualDuration = brent(f, interval[0], interval[1], 0.0001);
  } else if (!R_isnancpp(accrualDuration) && R_isnancpp(followupTime)) {
    // Lambda function
    auto f = [&](double t) {
      lr = lrpower(kMax, informationRates, criticalValues, futilityBounds,
                   allocationRatioPlanned, accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, stratumFraction,
                   lambda1, lambda2, gamma1, gamma2,
                   accrualDuration, t, fixedFollowup,
                   rho1, rho2, numSubintervals);
      return double(lr[0]) - (1-beta);
    };
    followupTime = brent(f, interval[0], interval[1], 0.0001);
  }

  // output the results
  lr = lrpower(kMax, informationRates, criticalValues, futilityBounds,
               allocationRatioPlanned, accrualTime, accrualIntensity,
               piecewiseSurvivalTime, stratumFraction,
               lambda1, lambda2, gamma1, gamma2,
               accrualDuration, followupTime, fixedFollowup,
               rho1, rho2, numSubintervals);
  return lr;
}
