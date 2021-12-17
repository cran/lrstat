#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Quantile function of truncated piecewise exponential distribution
//' @description Obtains the quantile of a piecewise expoenential distribution
//' given that it exceeds a specified lower bound.
//'
//' @param probability The scalar probability corresponding to the quantile.
//' @param piecewiseSurvivalTime A vector of starting times defining the time
//' pieces, must start with 0. Defaults to 0 for exponential distribution.
//' @param lambda A vector of hazard rates, one for each time interval.
//' @param lowerBound The left truncation time point for the survival time.
//' Defaults to zero for no truncation.
//'
//' @return The quantile x such that
//' P(X > x | X > lowerBound) = 1 - probability.
//'
//' @keywords internal
//'
//' @examples
//' qtpwexp1(probability = 0.3, piecewiseSurvivalTime = c(0, 6, 9, 15),
//'          lambda = c(0.025, 0.04, 0.015, 0.007))
//'
//' @export
// [[Rcpp::export]]
double qtpwexp1(const double probability = NA_REAL,
                const NumericVector& piecewiseSurvivalTime = 0,
                const NumericVector& lambda = NA_REAL,
                const double lowerBound = 0) {

  int j, j1, m;
  double q, v, v1;

  // cumulative hazard from lowerBound until the quantile
  v1 = -log(1 - probability);

  // identify the time interval containing the lowerBound
  m = piecewiseSurvivalTime.size();
  for (j=0; j<m; j++) {
    if (piecewiseSurvivalTime[j] > lowerBound) break;
  }
  j1 = (j==0 ? 0 : j-1); // to handle floating point precision

  if (j1 == m-1) { // in the last interval
    q = v1/lambda[j1] + lowerBound;
  } else {
    // accumulate the pieces on the cumulative hazard scale
    v = 0;
    for (j=j1; j<m-1; j++) {
      if (j==j1) {
        v += lambda[j]*(piecewiseSurvivalTime[j+1] - lowerBound);
      } else {
        v += lambda[j]*(piecewiseSurvivalTime[j+1] - piecewiseSurvivalTime[j]);
      }
      if (v >= v1) break;
    }

    if (j == m-1) { // in the last interval
      q = (v1 - v)/lambda[j] + piecewiseSurvivalTime[j];
    } else {
      q = piecewiseSurvivalTime[j+1] - (v - v1)/lambda[j];
    }
  }

  return q;
}


//' @title Quantile function of truncated piecewise exponential distribution
//' @description Obtains the quantile of a piecewise expoenential distribution
//' given that it exceeds a specified lower bound.
//'
//' @param probability The probabilities corresponding to the quantiles.
//' @param piecewiseSurvivalTime A vector of starting times defining the time
//' pieces, must start with 0. Defaults to 0 for exponential distribution.
//' @param lambda A vector of hazard rates, one for each time interval.
//' @param lowerBound The left truncation time point for the survival time.
//' Defaults to 0 for no truncation.
//'
//' @return A vector of quantiles x such that
//' P(X > x | X > lowerBound) = 1 - probability.
//'
//' @examples
//' qtpwexp(probability = c(0.3, 0.5), piecewiseSurvivalTime = c(0, 6, 9, 15),
//'         lambda = c(0.025, 0.04, 0.015, 0.007))
//'
//' @export
// [[Rcpp::export]]
NumericVector qtpwexp(const NumericVector& probability = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda = NA_REAL,
                      double lowerBound = 0) {

  int i, j, j1, k=probability.size(), m=piecewiseSurvivalTime.size();
  NumericVector v1(k);

  v1 = -log(1 - probability);

  // identify the time interval containing the lowerBound
  for (j=0; j<m; j++) {
    if (piecewiseSurvivalTime[j] > lowerBound) break;
  }
  j1 = (j==0 ? 0 : j-1);

  double v;
  NumericVector q(k);
  NumericVector t = piecewiseSurvivalTime;

  if (j1==m-1) { // in the last interval
    q = v1/lambda[j1] + lowerBound;
  } else {
    for (i=0; i<k; i++) {
      // accumulate the pieces on the cumulative hazard scale
      v = 0;
      for (j=j1; j<m-1; j++) {
        if (j==j1) {
          v += lambda[j]*(t[j+1] - lowerBound);
        } else {
          v += lambda[j]*(t[j+1] - t[j]);
        }
        if (v >= v1[i]) break;
      }

      if (j==m-1) { // in the last interval
        q[i] = (v1[i] - v)/lambda[j] + t[j];
      } else {
        q[i] = t[j+1] - (v - v1[i])/lambda[j];
      }
    }
  }

  return q;
}



//' @title Log-rank test simulation
//' @description Performs simulation for two-arm group sequential superiority
//' trials based on log-rank test.
//'
//' @param informationRates The information rates fixed prior to the trial.
//' Defaults to (1:kMax)/kMax if left unspecified.
//' @param kMax The maximum number of stages.
//' @param criticalValues Upper boundaries on the z-test statistic scale for
//' stopping for efficacy.
//' @param futilityBounds Lower boundaries on the z-test statistic scale
//' for stopping for futility at stages 1, ..., kMax-1. Defaults to
//' rep(-Inf, kMax-1) if left unspecified.
//' @param allocation1 Number of subjects in the active treatment group in
//' a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//' a randomization block. Defaults to 1 for equal randomization.
//' @param accrualTime Accrual time intervals, must start with 0, e.g.,
//' c(0, 3) breaks the time axis into 2 accrual intervals: [0, 3], (3, Inf).
//' Defaults to 0 for uniform accrual.
//' @param accrualIntensity A vector of accrual intensities, one for each
//' accrual time interval.
//' @param piecewiseSurvivalTime A vector that specifies the time intervals for
//' the piecewise exponential survival distribution, must start with 0, e.g.,
//' c(0, 6) breaks the time axis into 2 event intervals: [0, 6] and (6, Inf).
//' Defaults to 0 for exponential distribution.
//' @param stratumFraction A vector of stratum fractions.
//' Defaults to 1 for no stratification.
//' @param lambda1 A vector of hazard rates for the event for the
//' active treatment group, one for each analysis time interval, by stratum.
//' @param lambda2 A vector of hazard rates for the event for the
//' control group, one for each analysis time interval, by stratum.
//' @param gamma1 The hazard rate for exponential dropout or a vector of hazard
//' rates for piecewise exponential dropout for the active treatment group.
//' Defaults to 0 for no dropout.
//' @param gamma2 The hazard rate for exponential dropout or a vector of hazard
//' rates for piecewise exponential dropout for the control group.
//' Defaults to 0 for no dropout.
//' @param accrualDuration Duration of the enrollment period.
//' @param followupTime Follow-up time for the last enrolled subject.
//' @param fixedFollowup Whether a fixed follow-up design is used.
//' Defaults to 0 for variable follow-up.
//' @param rho1 First parameter of the Fleming-Harrington family of weighted
//' log-rank test. Defaults to 0 for conventional log-rank test.
//' @param rho2 Second parameter of the Fleming-Harrington family of weighted
//' log-rank test. Defaults to 0 for conventional log-rank test.
//' @param plannedEvents The planned cumulative total number of events at each
//' stage.
//' @param maxNumberOfIterations The number of simulation iterations.
//' Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//' to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//' The computer clock will be used if left unspecified,
//'
//' @return A list of S3 class lrsim with 3 components: sumstat is a list of
//' the operating characteristics of the design, sumdata is a data frame for
//' the summary data for each iteration, and rawdata is a data frame for
//' selected raw data if maxNumberOfRawDatasetsPerStage is a positive integer.
//'
//' @examples
//' sim = lrsim(kMax = 2, informationRates = c(0.5, 1),
//'             criticalValues = c(2.797, 1.977),
//'             accrualIntensity = 11,
//'             lambda1 = 0.018, lambda2 = 0.030,
//'             accrualDuration = 12,
//'             plannedEvents = c(75, 150),
//'             maxNumberOfIterations = 1000,
//'             maxNumberOfRawDatasetsPerStage = 1,
//'             seed = 314159)
//'
//' # summary statistics
//' sim
//'
//' # summary for each simulated data set
//' head(sim$sumdata)
//'
//' # raw data for selected replication
//' head(sim$rawdata)
//'
//' @export
// [[Rcpp::export]]
List lrsim(const int kMax = NA_INTEGER,
           NumericVector informationRates = NA_REAL,
           const NumericVector& criticalValues = NA_REAL,
           NumericVector futilityBounds = NA_REAL,
           const int allocation1 = 1,
           const int allocation2 = 1,
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
           const IntegerVector& plannedEvents = NA_INTEGER,
           const int maxNumberOfIterations = 1000,
           const int maxNumberOfRawDatasetsPerStage = 0,
           int seed = NA_INTEGER) {

  // b1 and b2 are the available slots for the two treatments in a block
  int i, iter, j, j1, j2, k, l, h, nsub;
  int nstrata = stratumFraction.size() ;
  double u, enrollt;

  IntegerVector b1(nstrata);
  IntegerVector b2(nstrata);

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


  NumericVector t = informationRates;

  NumericVector w(kMax);
  for (k=0; k<kMax; k++) {
    if (k==0) {
      w[k] = sqrt(t[k]);
    } else {
      w[k] = sqrt(t[k] - t[k-1]);
    }
  }

  // maximum number of subjects to enroll
  int m = accrualTime.size();
  double s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  int n = floor(s + 0.5);


  IntegerVector subjectId(n), stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  observationTime(n), timeUnderObservation(n), totalTime(n);

  LogicalVector event(n), dropoutEvent(n);


  IntegerVector repl(n), look(n), subj(n), strat(n), treat(n);

  NumericVector atime(n), survt(n), dropt(n), obst(n),
  ftime(n), totalt(n), evtime(n);

  LogicalVector status(n), drop(n);


  NumericVector cumStratumFraction = cumsum(stratumFraction);

  int nintervals = piecewiseSurvivalTime.size();
  NumericVector lam1(nintervals), lam2(nintervals);

  int nevents, nstages;
  IntegerVector accruals1(kMax), accruals2(kMax), totalAccruals(kMax),
  events1(kMax), events2(kMax), totalEvents(kMax);

  NumericVector analysisTime(kMax);

  int bigd = plannedEvents[kMax], bigd1;

  int d1, d2, dt, n1, n2, nt;

  // event and at-risk flags
  LogicalVector flag(bigd);

  double uscore1, vscore1, km1, w1;
  NumericVector uscore(kMax), vscore(kMax), lrstat(kMax), theta(kMax),
  zstat1(kMax), zstat(kMax);

  LogicalVector rejectPerStage(kMax), futilityPerStage(kMax);
  int stopStage;

  LogicalVector sub(n);

  // cache for the remaining number of raw data sets per stage to extract
  IntegerVector niter = rep(maxNumberOfRawDatasetsPerStage, kMax);

  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumberx(nrow1);
  IntegerVector stopStagex(nrow1);
  IntegerVector subjectIdx(nrow1);
  NumericVector arrivalTimex(nrow1);
  IntegerVector stratumx(nrow1);
  IntegerVector treatmentGroupx(nrow1);
  NumericVector survivalTimex(nrow1);
  NumericVector dropoutTimex(nrow1);
  NumericVector observationTimex(nrow1);
  NumericVector timeUnderObservationx(nrow1);
  LogicalVector eventx(nrow1);
  LogicalVector dropoutEventx(nrow1);

  int nrow2 = kMax*maxNumberOfIterations;

  IntegerVector iterationNumbery(nrow2);
  IntegerVector stageNumbery(nrow2);
  NumericVector analysisTimey(nrow2);
  IntegerVector accruals1y(nrow2);
  IntegerVector accruals2y(nrow2);
  IntegerVector totalAccrualsy(nrow2);
  IntegerVector events1y(nrow2);
  IntegerVector events2y(nrow2);
  IntegerVector totalEventsy(nrow2);
  NumericVector uscorey(nrow2);
  NumericVector vscorey(nrow2);
  NumericVector logRankStatisticy(nrow2);
  NumericVector hazardRatioEstimateLRy(nrow2);
  NumericVector zStatisticCHWy(nrow2);
  LogicalVector rejectPerStagey(nrow2);
  LogicalVector futilityPerStagey(nrow2);

  IntegerVector idx(kMax);

  DataFrame rawdata, sumdata;
  List sumstat, result;


  int nadded = 0;


  // set up random seed;
  if (seed==NA_INTEGER) {
    set_seed(std::time(0));
  } else {
    set_seed(seed);
  }


  for (iter=0; iter<maxNumberOfIterations; iter++) {
    b1.fill(allocation1);
    b2.fill(allocation2);

    enrollt = 0;
    for (i=0; i<n; i++) {
      subjectId[i] = i+1;


      // generate stratum information
      u = R::runif(0,1);
      for (j=0; j<nstrata; j++) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }


      // stratified block randomization;
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        b1[j] -= 1;
      } else {
        treatmentGroup[i] = 2;
        b2[j] -= 1;
      }

      // start a new block after depleting the current block;
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }


      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexp1(u, accrualTime, accrualIntensity, enrollt);
      arrivalTime[i] = enrollt;

      // generate survival time
      j1 = (stratum[i] - 1)*nintervals;
      j2 = stratum[i]*nintervals - 1;

      lam1 = lambda1[Range(j1,j2)];
      lam2 = lambda2[Range(j1,j2)];

      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexp1(u, piecewiseSurvivalTime, lam1);
      } else {
        survivalTime[i] = qtpwexp1(u, piecewiseSurvivalTime, lam2);
      }

      // generate dropout time
      u = R::runif(0,1);
      if (gamma1.size() == nintervals) {
        if (treatmentGroup[i]==1) {
          dropoutTime[i] = qtpwexp1(u, piecewiseSurvivalTime, gamma1);
        } else {
          dropoutTime[i] = qtpwexp1(u, piecewiseSurvivalTime, gamma2);
        }
      } else {
        if (treatmentGroup[i]==1) {
          dropoutTime[i] = qtpwexp1(u, 0, gamma1);
        } else {
          dropoutTime[i] = qtpwexp1(u, 0, gamma2);
        }
      }

      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime[i] < dropoutTime[i] &&
            survivalTime[i] < followupTime) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else if (dropoutTime[i] < survivalTime[i] &&
          dropoutTime[i] < followupTime) {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        } else {
          timeUnderObservation[i] = followupTime;
          event[i] = 0;
          dropoutEvent[i] = 0;
        }
      } else { // variable follow-up design
        if (survivalTime[i] < dropoutTime[i]) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        }
      }

      totalTime[i] = arrivalTime[i] + timeUnderObservation[i];

    }

    // find the analysis time for each stage
    nevents = sum(event);
    totalt = stl_sort(totalTime[event==1]);
    nstages = kMax;

    for (j=0; j<kMax; j++) {
      if (plannedEvents[j] >= nevents) {
        nstages = j+1;
        break;
      }
    }


    if (j==kMax) {
      totalEvents = clone(plannedEvents);
      for (k=0; k<nstages; k++) {
        analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
      }
    } else {
      for (k=0; k<nstages; k++) {
        if (k < nstages-1) {
          totalEvents[k] = plannedEvents[k];
          analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
        } else {
          totalEvents[k] = nevents;
          analysisTime[k] = totalt[nevents-1] + 1e-12;
        }
      }
    }



    // construct the log-rank test statistic at each stage
    stopStage = nstages;
    for (k=0; k<nstages; k++) {

      // subset data included for the kth analysis
      sub = (arrivalTime < analysisTime[k]);
      nsub = sum(sub);

      repl = rep(iter+1,nsub);
      look = rep(k+1, nsub);
      subj = subjectId[sub];
      atime = arrivalTime[sub];
      strat = stratum[sub];
      treat = treatmentGroup[sub];
      survt = survivalTime[sub];
      dropt = dropoutTime[sub];
      obst = rep(analysisTime[k], nsub);
      ftime = timeUnderObservation[sub];
      status = event[sub];
      drop = dropoutEvent[sub];

      for (i=0; i<nsub; i++) {
        if (ftime[i] > analysisTime[k] - atime[i]) {
          ftime[i] = analysisTime[k] - atime[i];
          status[i] = 0;
          drop[i] = 0;
        }
      }


      // number of accrued patients and total number of events
      accruals1[k] = sum(treat==1);
      accruals2[k] = sum(treat==2);
      totalAccruals[k] = accruals1[k] + accruals2[k];

      events1[k] = sum((status==1)*(treat==1));
      events2[k] = sum((status==1)*(treat==2));
      bigd1 = totalEvents[k];

      // identify the event times
      evtime = stl_sort(ftime[status==1]);

      // count number of patients with event and number of patients at risk
      // by stratum and event time point
      uscore[k] = 0;
      vscore[k] = 0;
      for (h=1; h<=nstrata; h++) {
        // log-rank statistic in the stratum
        uscore1 = 0;
        vscore1 = 0;
        km1 = 1;

        for (j=0; j<bigd1; j++) {
          // Fleming-Harrington weight
          w1 = pow(km1, rho1)*pow(1-km1, rho2);

          // number of events at the time point
          flag = (strat==h) * (ftime==evtime[j]) * (status==1);
          d1 = sum( flag * (treat==1) );
          d2 = sum( flag * (treat==2) );
          dt = d1 + d2;

          // number of patients at risk at the time point
          flag = (strat==h) * (ftime>=evtime[j]);
          n1 = sum( flag * (treat==1) );
          n2 = sum( flag * (treat==2) );
          nt = n1 + n2;

          // accumulate the mean and variance of the score statistic
          if (nt > 1) {
            uscore1 += w1*(d1 - n1*dt/(nt+0.0));
            vscore1 += w1*w1*dt*(nt-dt)*n1*n2/(nt*nt*(nt-1)+0.0);
          }

          // update the KM estimate for the pooled sample
          km1 = km1*(1-dt/(nt+0.0));
        }
        // sum across strata
        uscore[k] += uscore1;
        vscore[k] += vscore1;
      }

      // log-rank z statistic and associated hazard ratio estimate
      lrstat[k] = uscore[k]/sqrt(vscore[k]);
      theta[k] = exp(uscore[k]/vscore[k]);

      // stage-specific z statistic assuming lower hazard ratio is better
      if (k==0) {
        zstat1[k] = -lrstat[k];
      } else {
        zstat1[k] = -(uscore[k] - uscore[k-1])/sqrt(vscore[k] - vscore[k-1]);
      }

      // CHW z statistic to account for information time different from planned
      if (k==0) {
        zstat[k] = zstat1[k];
      } else {
        zstat[k] = (sqrt(t[k-1])*zstat[k-1] + w[k]*zstat1[k]) / sqrt(t[k]);
      }


      // compare to the critical values to make decisions
      rejectPerStage[k] = 0;
      futilityPerStage[k] = 0;
      if (zstat[k] > criticalValues[k]) {
        rejectPerStage[k] = 1;
      } else if ((k < nstages-1 && zstat[k] < futilityBounds[k])
                   || (k == nstages-1))  {
        futilityPerStage[k] = 1;
      }

      if (rejectPerStage[k]==1 || futilityPerStage[k]==1) {

        // add raw data to output
        if (niter[k] > 0) {
          if (nadded==0) { // initialize
            iterationNumberx = clone(repl);
            stopStagex = clone(look);
            subjectIdx = clone(subj);
            arrivalTimex = clone(atime);
            stratumx = clone(strat);
            treatmentGroupx = clone(treat);
            survivalTimex = clone(survt);
            dropoutTimex = clone(dropt);
            observationTimex = clone(obst);
            timeUnderObservationx = clone(ftime);
            eventx = clone(status);
            dropoutEventx = clone(drop);
          } else { // append
            for (l=0; l<nsub; l++) {
              iterationNumberx.push_back(repl[l]);
              stopStagex.push_back(look[l]);
              subjectIdx.push_back(subj[l]);
              arrivalTimex.push_back(atime[l]);
              stratumx.push_back(strat[l]);
              treatmentGroupx.push_back(treat[l]);
              survivalTimex.push_back(survt[l]);
              dropoutTimex.push_back(dropt[l]);
              observationTimex.push_back(obst[l]);
              timeUnderObservationx.push_back(ftime[l]);
              eventx.push_back(status[l]);
              dropoutEventx.push_back(drop[l]);
            }
          }
          // update the number of stage k dataset to extract
          niter[k] = niter[k]-1;
          nadded++;
        }

        stopStage = k+1;
        break;
      }


    }


    // add summary data to output
    if (iter==0) {
      idx = Range(0, stopStage-1);
      iterationNumbery = rep(iter+1, stopStage);
      stageNumbery = seq_len(stopStage);
      analysisTimey = analysisTime[idx];
      accruals1y = accruals1[idx];
      accruals2y = accruals2[idx];
      totalAccrualsy = totalAccruals[idx];
      events1y = events1[idx];
      events2y = events2[idx];
      totalEventsy = totalEvents[idx];
      uscorey = uscore[idx];
      vscorey = vscore[idx];
      logRankStatisticy = lrstat[idx];
      hazardRatioEstimateLRy = theta[idx];
      zStatisticCHWy = zstat[idx];
      rejectPerStagey = rejectPerStage[idx];
      futilityPerStagey = futilityPerStage[idx];
    } else {
      for (l=0; l<stopStage; l++) {
        iterationNumbery.push_back(iter+1);
        stageNumbery.push_back(l+1);
        analysisTimey.push_back(analysisTime[l]);
        accruals1y.push_back(accruals1[l]);
        accruals2y.push_back(accruals2[l]);
        totalAccrualsy.push_back(totalAccruals[l]);
        events1y.push_back(events1[l]);
        events2y.push_back(events2[l]);
        totalEventsy.push_back(totalEvents[l]);
        uscorey.push_back(uscore[l]);
        vscorey.push_back(vscore[l]);
        logRankStatisticy.push_back(lrstat[l]);
        hazardRatioEstimateLRy.push_back(theta[l]);
        zStatisticCHWy.push_back(zstat[l]);
        rejectPerStagey.push_back(rejectPerStage[l]);
        futilityPerStagey.push_back(futilityPerStage[l]);
      }
    }

  }


  // simulation results on power and expected sample size

  NumericVector pRejectPerStage(kMax), pFutilityPerStage(kMax),
  nEventsPerStage(kMax), nSubjectsPerStage(kMax),
  analysisTimePerStage(kMax);


  // number of observations in the summary dataset
  int nrow3 = stageNumbery.size();

  // number of iterations for each stage
  IntegerVector count(kMax);

  for (i=0; i<nrow3; i++) {
    k = stageNumbery[i] - 1;
    count[k] += 1;
    pRejectPerStage[k] += rejectPerStagey[i];
    pFutilityPerStage[k] += futilityPerStagey[i];
    nEventsPerStage[k] += totalEventsy[i];
    nSubjectsPerStage[k] += totalAccrualsy[i];
    analysisTimePerStage[k] += analysisTimey[i];
  }

  for (k=0; k<kMax; k++) {
    if (count[k] > 0) {
      pRejectPerStage[k] /= maxNumberOfIterations;
      pFutilityPerStage[k] /= maxNumberOfIterations;
      nEventsPerStage[k] /= count[k];
      nSubjectsPerStage[k] /= count[k];
      analysisTimePerStage[k] /= count[k];
    }
  }


  double pOverallReject = sum(pRejectPerStage);

  double expectedNumberOfEvents=0, expectedNumberOfSubjects=0,
  expectedStudyDuration=0;

  iter = 1;
  k = 1;
  for (i=0; i<nrow3; i++) {
    if (iterationNumbery[i] > iter) {
      // accumulate at each iteration
      expectedNumberOfEvents += totalEventsy[i-1];
      expectedNumberOfSubjects += totalAccrualsy[i-1];
      expectedStudyDuration += analysisTimey[i-1];
      // advance the iteration and reset the stage
      iter++;
      k = 1;
    } else if (stageNumbery[i] > k) {
      // march toward the stop stage of the iteration
      k++;
    }
  }
  // add the information for the last iteration
  expectedNumberOfEvents += totalEventsy[i-1];
  expectedNumberOfSubjects += totalAccrualsy[i-1];
  expectedStudyDuration += analysisTimey[i-1];

  // take the mean
  expectedNumberOfEvents /= maxNumberOfIterations;
  expectedNumberOfSubjects /= maxNumberOfIterations;
  expectedStudyDuration /= maxNumberOfIterations;


  sumstat = List::create(_["overallReject"]=pOverallReject,
                         _["rejectPerStage"]=pRejectPerStage,
                         _["futilityPerStage"]=pFutilityPerStage,
                         _["eventsPerStage"]=nEventsPerStage,
                         _["numberOfSubjects"]=nSubjectsPerStage,
                         _["analysisTime"]=analysisTimePerStage,
                         _["expectedNumberOfEvents"]=expectedNumberOfEvents,
                         _["expectedNumberOfSubjects"]=expectedNumberOfSubjects,
                         _["expectedStudyDuration"]=expectedStudyDuration);




  // simulation datasets

  sumdata = DataFrame::create(_["iterationNumber"]=iterationNumbery,
                              _["stageNumber"]=stageNumbery,
                              _["analysisTime"]=analysisTimey,
                              _["accruals1"]=accruals1y,
                              _["accruals2"]=accruals2y,
                              _["totalAccruals"]=totalAccrualsy,
                              _["events1"]=events1y,
                              _["events2"]=events2y,
                              _["totalEvents"]=totalEventsy,
                              _["uscore"]=uscorey,
                              _["vscore"]=vscorey,
                              _["logRankStatistic"]=logRankStatisticy,
                              _["hazardRatioEstimateLR"]=hazardRatioEstimateLRy,
                              _["zStatisticCHW"]=zStatisticCHWy,
                              _["rejectPerStage"]=rejectPerStagey,
                              _["futilityPerStage"]=futilityPerStagey);



  if (maxNumberOfRawDatasetsPerStage > 0) {
    rawdata = DataFrame::create(_["iterationNumber"]=iterationNumberx,
                                _["stopStage"]=stopStagex,
                                _["subjectId"]=subjectIdx,
                                _["arrivalTime"]=arrivalTimex,
                                _["stratum"]=stratumx,
                                _["treatmentGroup"]=treatmentGroupx,
                                _["survivalTime"]=survivalTimex,
                                _["dropoutTime"]=dropoutTimex,
                                _["observationTime"]=observationTimex,
                                _["timeUnderObservation"]=timeUnderObservationx,
                                _["event"]=eventx,
                                _["dropoutEvent"]=dropoutEventx);

    result = List::create(_["sumstat"]=sumstat,
                          _["sumdata"]=sumdata,
                          _["rawdata"]=rawdata);
  } else {
    result = List::create(_["sumstat"]=sumstat,
                          _["sumdata"]=sumdata);
  }

  result.attr("class") = "lrsim";


  return result;

}



