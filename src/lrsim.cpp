#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;



//' @title Log-rank test simulation
//' @description Performs simulation for two-arm group sequential superiority
//' trials based on log-rank test.
//'
//' @inheritParams param_kMax
//' @param informationTime Information time in terms of variance of
//'   weighted log-rank test score statistic. Same as informationRates
//'   in terms of number of events for unweighted log-rank test. Use caltime
//'   and lrstat to find the information time for weighted log-rank tests.
//'   Fixed prior to the trial. Defaults to \code{(1:kMax) / kMax} if
//'   left unspecified.
//' @inheritParams param_criticalValues
//' @inheritParams param_futilityBounds
//' @param allocation1 Number of subjects in the active treatment group in
//' a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//' a randomization block. Defaults to 1 for equal randomization.
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
//' @param plannedEvents The planned cumulative total number of events at each
//' stage.
//' @param maxNumberOfIterations The number of simulation iterations.
//' Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//' to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//' The computer clock will be used if left unspecified,
//'
//' @return A list of S3 class lrsim with 3 components: overview is a list of
//' the operating characteristics of the design, sumdata is a data frame for
//' the summary data for each iteration, and rawdata is a data frame for
//' selected raw data if maxNumberOfRawDatasetsPerStage is a positive integer.
//'
//' @examples
//' sim = lrsim(kMax = 2, informationTime = c(0.5, 1),
//'             criticalValues = c(2.797, 1.977),
//'             accrualIntensity = 11,
//'             lambda1 = 0.018, lambda2 = 0.030,
//'             accrualDuration = 12,
//'             plannedEvents = c(60, 120),
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
           NumericVector informationTime = NA_REAL,
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

  int i, iter, j, j1, j2, k, h, nsub;
  int nstrata = stratumFraction.size() ;
  int nintervals = piecewiseSurvivalTime.size();
  double u, enrollt;

  // b1 and b2 are the available slots for the two treatments in a block
  IntegerVector b1(nstrata);
  IntegerVector b2(nstrata);

  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  } else if (kMax <= 0) {
    stop("kMax must be a positive integer");
  }

  // set default parameter values
  if (is_false(any(is_na(informationTime)))) {
    if (informationTime.size() != kMax) {
      stop("Invalid length for informationTime");
    } else if (informationTime[0] <= 0) {
      stop("Elements of informationTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationTime) <= 0))) {
      stop("Elements of informationTime must be increasing");
    } else if (informationTime[kMax-1] != 1) {
      stop("informationTime must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationTime = as<NumericVector>(tem)/(kMax+0.0);
  }


  if (is_true(any(is_na(criticalValues)))) {
    stop("criticalValues must be provided");
  } else if (criticalValues.size() != kMax) {
    stop("Invalid length for criticalValues");
  }


  if (kMax > 1) {
    if (is_true(any(is_na(futilityBounds)))) {
      futilityBounds = rep(R_NegInf, kMax-1);
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }
  }

  if (allocation1 <= 0) {
    stop("allocation1 must be positive");
  }

  if (allocation2 <= 0) {
    stop("allocation2 must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  } else if (kMax > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  } else if (kMax > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  } else if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }


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

  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  } else if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    } else if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }


  if (plannedEvents[0] <= 0) {
    stop("Elements of plannedEvents must be positive");
  } else if (plannedEvents.size() > 1 &&
    is_true(any(diff(plannedEvents) <= 0))) {
    stop("Elements of plannedEvents must be increasing");
  }

  if (maxNumberOfIterations <= 0) {
    stop("maxNumberOfIterations must be positive");
  }

  if (maxNumberOfRawDatasetsPerStage  < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be non-negative");
  }


  NumericVector t = informationTime;

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


  IntegerVector stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  observationTime(n), timeUnderObservation(n), totalTime(n), totalt(n);

  LogicalVector event(n), dropoutEvent(n);

  NumericVector cumStratumFraction = cumsum(stratumFraction);

  NumericVector lam1(nintervals), lam2(nintervals);

  int nevents, nstages;
  NumericVector analysisTime(kMax);

  IntegerVector accruals1(kMax), accruals2(kMax), totalAccruals(kMax),
  events1(kMax), events2(kMax), totalEvents(kMax);

  NumericVector timeUnderObservationSorted(n);
  IntegerVector sortedIndex(n), stratumSorted(n), treatmentGroupSorted(n);
  LogicalVector eventSorted(n);

  double uscore1, vscore1;
  NumericVector km1(nstrata), w1(nstrata);

  NumericVector uscore(kMax), vscore(kMax), lrstat(kMax),
  zstat1(kMax), zstat(kMax);

  LogicalVector rejectPerStage(kMax), futilityPerStage(kMax);
  int stopStage;

  LogicalVector sub(n);

  IntegerVector n1(nstrata), n2(nstrata), nt(nstrata);

  // cache for the number of raw data sets per stage to extract
  IntegerVector niter(kMax);

  int nrow1 = std::min(n*kMax*maxNumberOfRawDatasetsPerStage,
                       n*maxNumberOfIterations);

  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stopStagex = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector subjectIdx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector arrivalTimex = NumericVector(nrow1, NA_REAL);
  IntegerVector stratumx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector treatmentGroupx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector survivalTimex = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTimex = NumericVector(nrow1, NA_REAL);
  NumericVector observationTimex = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservationx = NumericVector(nrow1, NA_REAL);
  LogicalVector eventx = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEventx = LogicalVector(nrow1, NA_LOGICAL);

  int nrow2 = kMax*maxNumberOfIterations;

  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector stageNumbery = IntegerVector(nrow2, NA_INTEGER);
  NumericVector analysisTimey = NumericVector(nrow2, NA_REAL);
  IntegerVector accruals1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalAccrualsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalEventsy = IntegerVector(nrow2, NA_INTEGER);
  NumericVector uscorey = NumericVector(nrow2, NA_REAL);
  NumericVector vscorey = NumericVector(nrow2, NA_REAL);
  NumericVector logRankStatisticy = NumericVector(nrow2, NA_REAL);
  NumericVector zStatisticCHWy = NumericVector(nrow2, NA_REAL);
  LogicalVector rejectPerStagey = LogicalVector(nrow2, NA_LOGICAL);
  LogicalVector futilityPerStagey = LogicalVector(nrow2, NA_LOGICAL);



  int index1=0, index2=0;
  double time;

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
      enrollt = qtpwexp(u, accrualTime, accrualIntensity, enrollt);
      arrivalTime[i] = enrollt;

      // generate survival time
      j1 = j*nintervals;
      j2 = j1 + nintervals - 1;

      lam1 = lambda1[Range(j1,j2)];
      lam2 = lambda2[Range(j1,j2)];

      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexp(u, piecewiseSurvivalTime, lam1, 0);
      } else {
        survivalTime[i] = qtpwexp(u, piecewiseSurvivalTime, lam2, 0);
      }

      // generate dropout time
      u = R::runif(0,1);
      if (gamma1.size() == nintervals) {
        if (treatmentGroup[i]==1) {
          dropoutTime[i] = qtpwexp(u, piecewiseSurvivalTime, gamma1, 0);
        } else {
          dropoutTime[i] = qtpwexp(u, piecewiseSurvivalTime, gamma2, 0);
        }
      } else {
        if (treatmentGroup[i]==1) {
          dropoutTime[i] = qtpwexp(u, 0, gamma1, 0);
        } else {
          dropoutTime[i] = qtpwexp(u, 0, gamma2, 0);
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
    totalt = stl_sort(totalTime[event]);
    nstages = kMax;

    for (j=0; j<kMax; j++) {
      if (plannedEvents[j] >= nevents) {
        nstages = j+1;
        break;
      }
    }


    if (j==kMax) { // total number of events exceeds planned
      for (k=0; k<nstages; k++) {
        analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
      }
    } else {
      for (k=0; k<nstages; k++) {
        if (k < nstages-1) {
          analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
        } else {
          analysisTime[k] = totalt[nevents-1] + 1e-12;
        }
      }
    }

    // construct the log-rank test statistic at each stage
    stopStage = nstages;
    for (k=0; k<nstages; k++) {
      time = analysisTime[k];

      n1.fill(0);  // number of subjects in each stratum by treatment
      n2.fill(0);
      events1[k] = 0;
      events2[k] = 0;

      for (i=0; i<n; i++) {
        h = stratum[i]-1;
        observationTime[i] = time;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation[i] = time - arrivalTime[i];
          event[i] = 0;
          dropoutEvent[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            n1[h]++;
          } else if (treatmentGroup[i]==2) {
            n2[h]++;
          }

          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime[i] < time &&
                survivalTime[i] < dropoutTime[i] &&
                survivalTime[i] < followupTime) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] < time &&
              dropoutTime[i] < survivalTime[i] &&
              dropoutTime[i] < followupTime) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime[i] &&
              followupTime < dropoutTime[i]) {
              timeUnderObservation[i] = followupTime;
              event[i] = 0;
              dropoutEvent[i] = 0;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime[i] < time &&
                survivalTime[i] < dropoutTime[i]) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] < time &&
              dropoutTime[i] < survivalTime[i]) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event[i]) events1[k]++;
          if (treatmentGroup[i]==2 && event[i]) events2[k]++;
        }
      }

      // number of accrued patients and total number of events
      accruals1[k] = 0;
      accruals2[k] = 0;
      for (h=0; h<nstrata; h++) {
        accruals1[k] += n1[h];
        accruals2[k] += n2[h];
      }
      totalAccruals[k] = accruals1[k] + accruals2[k];

      totalEvents[k] = events1[k] + events2[k];

      // order the data by time under observation
      timeUnderObservationSorted = stl_sort(timeUnderObservation);
      sortedIndex = match(timeUnderObservationSorted, timeUnderObservation);
      sortedIndex = sortedIndex - 1;
      eventSorted = event[sortedIndex];
      stratumSorted = stratum[sortedIndex];
      treatmentGroupSorted = treatmentGroup[sortedIndex];
      sub = (timeUnderObservationSorted > 0);
      eventSorted = eventSorted[sub];
      stratumSorted = stratumSorted[sub];
      treatmentGroupSorted = treatmentGroupSorted[sub];
      nsub = eventSorted.size();

      // calculate the stratified log-rank test
      uscore1 = 0;
      vscore1 = 0;
      km1.fill(1);  // km(t-) estimate by stratum
      for (i=0; i<nsub; i++) {
        h = stratumSorted[i] - 1;
        nt[h] = n1[h] + n2[h];

        if (eventSorted[i]) { // at most 1 event can occur at a given time
          w1[h] = pow(km1[h], rho1)*pow(1-km1[h], rho2);
          if (nt[h] > 0) {
            uscore1 += w1[h]*((treatmentGroupSorted[i]==1) - n1[h]/(nt[h]+0.0));
            vscore1 += w1[h]*w1[h]*n1[h]*n2[h]/(nt[h]*nt[h]+0.0);
          }
          km1[h] *= (1-1/(nt[h]+0.0)); // update km estimate
        }

				// reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          n1[h]--;
        } else {
          n2[h]--;
        }
      }

      uscore[k] = uscore1;
      vscore[k] = vscore1;

      // log-rank z statistic
      lrstat[k] = uscore[k]/sqrt(vscore[k]);

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
        if (niter[k] < maxNumberOfRawDatasetsPerStage) {
          for (i=0; i<n; i++) {
            iterationNumberx[index1] = iter+1;
            stopStagex[index1] = k+1;
            subjectIdx[index1] = i+1;
            arrivalTimex[index1] = arrivalTime[i];
            stratumx[index1] = stratum[i];
            treatmentGroupx[index1] = treatmentGroup[i];
            survivalTimex[index1] = survivalTime[i];
            dropoutTimex[index1] = dropoutTime[i];
            observationTimex[index1] = observationTime[i];
            timeUnderObservationx[index1] = timeUnderObservation[i];
            eventx[index1] = event[i];
            dropoutEventx[index1] = dropoutEvent[i];
            index1++;
          }

          // update the number of stage k dataset to extract
          niter[k]++;
        }

        stopStage = k+1;
        break;
      }

    }

    // add summary data to output
    for (k=0; k<stopStage; k++) {
      iterationNumbery[index2] = iter+1;
      stageNumbery[index2] = k+1;
      analysisTimey[index2] = analysisTime[k];
      accruals1y[index2] = accruals1[k];
      accruals2y[index2] = accruals2[k];
      totalAccrualsy[index2] = totalAccruals[k];
      events1y[index2] = events1[k];
      events2y[index2] = events2[k];
      totalEventsy[index2] = totalEvents[k];
      uscorey[index2] = uscore[k];
      vscorey[index2] = vscore[k];
      logRankStatisticy[index2] = lrstat[k];
      zStatisticCHWy[index2] = zstat[k];
      rejectPerStagey[index2] = rejectPerStage[k];
      futilityPerStagey[index2] = futilityPerStage[k];
      index2++;
    }


  }

  // only keep nonmissing records
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  stageNumbery = stageNumbery[sub2];
  analysisTimey = analysisTimey[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  totalEventsy = totalEventsy[sub2];
  uscorey = uscorey[sub2];
  vscorey = vscorey[sub2];
  logRankStatisticy = logRankStatisticy[sub2];
  zStatisticCHWy = zStatisticCHWy[sub2];
  rejectPerStagey = rejectPerStagey[sub2];
  futilityPerStagey = futilityPerStagey[sub2];



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
    count[k]++;
    pRejectPerStage[k] += rejectPerStagey[i];
    pFutilityPerStage[k] += futilityPerStagey[i];
    nEventsPerStage[k] += totalEventsy[i];
    nSubjectsPerStage[k] += totalAccrualsy[i];
    analysisTimePerStage[k] += analysisTimey[i];
  }

  // nEventsPerStage, nSubjectsPerStage, analysisTimePerStage are conditional
  // expectations given that the trial stops at the current or future stages
  for (k=0; k<kMax; k++) {
    if (count[k] > 0) {
      pRejectPerStage[k] /= maxNumberOfIterations;
      pFutilityPerStage[k] /= maxNumberOfIterations;
      nEventsPerStage[k] /= count[k];
      nSubjectsPerStage[k] /= count[k];
      analysisTimePerStage[k] /= count[k];
    }
  }

  NumericVector cpu = cumsum(pRejectPerStage);
  NumericVector cpl = cumsum(pFutilityPerStage);

  double pOverallReject = sum(pRejectPerStage);

  double expectedNumberOfEvents=0, expectedNumberOfSubjects=0,
    expectedStudyDuration=0;

  iter = 1;
  for (i=0; i<nrow3; i++) {
    if (iterationNumbery[i] > iter) {
      // accumulate at each iteration
      expectedNumberOfEvents += totalEventsy[i-1];
      expectedNumberOfSubjects += totalAccrualsy[i-1];
      expectedStudyDuration += analysisTimey[i-1];
      // advance the iteration
      iter++;
    }
  }
  // add the information for the last iteration
  expectedNumberOfEvents += totalEventsy[i-1];
  expectedNumberOfSubjects += totalAccrualsy[i-1];
  expectedStudyDuration += analysisTimey[i-1];

  // expectedNumberOfEvents, expectedNumberOfSubjects, expectedStudyDuration
  // are at the averages when the trial stops
  expectedNumberOfEvents /= maxNumberOfIterations;
  expectedNumberOfSubjects /= maxNumberOfIterations;
  expectedStudyDuration /= maxNumberOfIterations;

  List overview = List::create(
    _["rejectPerStage"]=pRejectPerStage,
    _["futilityPerStage"]=pFutilityPerStage,
    _["cumulativeRejection"]=cpu,
    _["cumulativeFutility"]=cpl,
    _["numberOfEvents"]=nEventsPerStage,
    _["numberOfSubjects"]=nSubjectsPerStage,
    _["analysisTime"]=analysisTimePerStage,
    _["overallReject"]=pOverallReject,
    _["expectedNumberOfEvents"]=expectedNumberOfEvents,
    _["expectedNumberOfSubjects"]=expectedNumberOfSubjects,
    _["expectedStudyDuration"]=expectedStudyDuration);



  // simulation datasets
  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"]=iterationNumbery,
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
    _["zStatisticCHW"]=zStatisticCHWy,
    _["rejectPerStage"]=rejectPerStagey,
    _["futilityPerStage"]=futilityPerStagey);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stopStagex = stopStagex[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTimex = survivalTimex[sub1];
    dropoutTimex = dropoutTimex[sub1];
    observationTimex = observationTimex[sub1];
    timeUnderObservationx = timeUnderObservationx[sub1];
    eventx = eventx[sub1];
    dropoutEventx = dropoutEventx[sub1];

    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"]=iterationNumberx,
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

    result = List::create(_["overview"]=overview,
                          _["sumdata"]=sumdata,
                          _["rawdata"]=rawdata);
  } else {
    result = List::create(_["overview"]=overview,
                          _["sumdata"]=sumdata);
  }

  result.attr("class") = "lrsim";


  return result;
}



