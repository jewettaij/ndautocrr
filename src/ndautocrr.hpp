///   @file ndautocrr.hpp
///   @brief  Calculate the auto-correlation function from multiple time series.
///   @date 2007-4-13

#ifndef _NDAUTOCRR_HPP
#define _NDAUTOCRR_HPP
#include <vector>
#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>
#include "err.h"
#include "inner_product.h"
using namespace std;


/// @brief This class is used to read data from multiple independent data sets
///        (arrays storing time series) and calculate the correlation function,
///        averaged over all of these data sets.  (The data contained
///        in each set is assumed to be independent of the others.)

template<typename Scalar>

class NdAutocrr {

  bool is_periodic;
  size_t L;
  Scalar threshold; //used for deciding when to cut the autocorrelation function
  Scalar persistence_length_threshold; //used for calculating persistence length
  bool subtract_ave;
  bool report_rms;

public:

  /// vC[j] stores the (discretized) correlation function (C(j))
  /// (The caller must read in the data and invoke Finalize() beforehand.)
  vector<double> vC;

  /// vCrms[j] stores the fluctuations around the mean value of vC[j]
  /// (The caller must read in the data and invoke Finalize() beforehand.)
  /// This quantity is only calculated when report_rms == true
  /// (and is probably not useful).
  vector<double> vCrms;

  /// vNumSamples[j]=number of terms averaged together
  /// to compute vC[j] (considering all data sets)
  /// (The caller must read in the data beforehand.)
  vector<size_t> vNumSamples; 

  NdAutocrr(Scalar _threshold=-1.01, //!< the value below which vC[j]/vc[0] must fall before it is discarded
            size_t _L=0, //!< _L+1 = the requested size of vC
            bool _is_periodic = false, //!< wrap i+j back into [0,N) when calculating x(i+j)?
            bool _subtract_ave=true, //!< Compute <(x(i)-<x>)*(x(i+j)-<x>)> OR <x(i)*x(i+j)> ?
            bool _report_rms=false //!< Calculate the rms values of (x(i)-<x>)*(x(i+j)-<x>) ?  (probably not useful)
            ):
    threshold(_threshold),
    L(_L),
    is_periodic(_is_periodic),
    subtract_ave(_subtract_ave),
    report_rms(_report_rms)
  {
    persistence_length_threshold = threshold;
    if (L > 0) {
      // If the user manually specified a domainwidth using the \"-L\" argument
      // then do not use thresholding to make the computation faster
      // (but do use it to calculate the persistence length).
      threshold = -1.01;
      vC.resize(L+1);
      if (report_rms)
        vCrms.resize(L+1);
      vNumSamples.resize(L+1);
    }
    else if (threshold <= -1.0) {
      threshold = 1.0 / M_E; //default threshold is 1/e
      persistence_length_threshold = threshold;
    }
  }


  /// @brief  Return the size of the vC[j] array (or at least the portion
  ///         of which we care about).
  size_t size() const { assert(L+1 == vC.size()); return L; }


  /// @brief  Accumulate the sums used to calculate the average (vC[j])
  size_t
  Accumulate(const vector<vector<vector<Scalar> > > &vvvX_nid, //!< series of data points (each is a vector of dimension d)
             ostream *pReportProgress = nullptr)  //!< print progress to the user?
  {
    for (int n=0; n < vvvX_nid.size(); n++)
      AccumulateSingle(vvvX_nid[n], pReportProgress);
    return L;
  }


  /// @brief Accumulate the sums used to calculate the average (vC[j])
  ///        considering only data from a single data set.
  size_t
  AccumulateSingle(const vector<vector<Scalar> > &vvX_id, //!< series of data points (each is a vector of dimension d)
                   ostream *pReportProgress = nullptr)  //!< print progress to the user?
  {
    vector<vector<Scalar> > _vvX_id = vvX_id; //make a local copy of the data set
    size_t L_single = _vvX_id.size();
    ChooseL(L_single);

    vector<Scalar> x_ave;

    int D = 0;
    for (size_t i=0; i < _vvX_id.size(); i++) {
      if (i == 0) 
        D = _vvX_id[i].size();
      else if (D != _vvX_id[i].size())
        throw InputErr("Error: Inconsistent number of entries on each line.\n");
    }


    // Should we subtract the average value when calculating autocorrelation?
    if (subtract_ave) {

      // If so, then ...
      x_ave.resize(D, 0.0);
      for (size_t i=0; i < _vvX_id.size(); i++)
        for (size_t d=0; d < D; d++)
          x_ave[d] += _vvX_id[i][d];
      for (size_t d=0; d < D; d++)
        x_ave[d] /= _vvX_id.size();

      // subtract the average value
      for (size_t i=0; i < _vvX_id.size(); ++i)
        for (size_t d=0; d < D; d++)
        _vvX_id[i][d] -= x_ave[d];
    }


    if (is_periodic)
    {

      #pragma omp parallel
      {
        #pragma omp for collapse(1)
        for (size_t j=0; j <= L; ++j)
        {
          if (pReportProgress)
            *pReportProgress << "#    processing separation " << j << endl;
          for (size_t i=0; i < _vvX_id.size(); ++i)
          {
            size_t iplusj = i+j;
            if (iplusj >= _vvX_id.size()) 
              iplusj -= _vvX_id.size();
            assert((0 <= iplusj) && (iplusj < _vvX_id.size()));

            Scalar C = inner_product(_vvX_id[i], _vvX_id[iplusj]);
            vC[j] += C;

            if (report_rms)
              vCrms[j] += C*C;
          }

          vNumSamples[j] += _vvX_id.size();

          // Check for threshold violations.
          // If the covariance function is too low, then quit
          if (vC[j] < threshold * vC[0])
          {
            // This will get us out of the loop, and also
            // truncate the domain of the covariance function
            L = j-1;
          }
        } //for (size_t j=0; j <= L; ++j)
      } //#pragma omp parallel

    } //if (is_periodic)
    else
    {

      size_t jmax = _vvX_id.size();
      if (jmax > L)
        jmax = L;

      #pragma omp parallel
      {
        #pragma omp for collapse(1)
        for (size_t j=0; j <= jmax; ++j)
        {
          if (j > jmax)
            continue;
          if (pReportProgress)
            *pReportProgress << "#    processing separation " << j << endl;
          for (size_t i=0; i < _vvX_id.size()-j; ++i)
          {
            Scalar C = inner_product(_vvX_id[i], _vvX_id[i+j]);
            vC[j] += C;

            if (report_rms)
              vCrms[j] += C*C;
          }

          vNumSamples[j] += _vvX_id.size()-j;

          // Check for threshold violations.
          // If the covariance function is too low, then quit
          if (vC[j] < threshold * vC[0]) {
            #pragma omp critical
            {
              if (j < jmax) {
                L = j;       //This will truncate the correlation function.
                jmax=j;      //This will break us out of the loop.
              }
            }
          }
        } //for (size_t j=0; j <= jmax; ++j)
      } //#pragma omp parallel

    } //else clause for "if (is_periodic)"

    if (vC.size() <= L) // if we reduced L, truncate the correlation function
      vC.resize(L+1);
    return L;
  } //size_t AccumulateSingle()




  /// @brief Invoke this function after reading all the data sets.
  ///        This function (for each separation length, j) divides vC and vCrms
  ///        by the total number of samples (accross all data sets, for that j).
  void
  Finalize() {

    assert(L+1 <= vC.size());

    if (L+1 < vC.size())
      vC.resize(L+1);
    if (L+1 < vCrms.size())
      vCrms.resize(L+1);

    for (size_t j=0; j < L+1; ++j) {
      if (vNumSamples[j] > 0) {
        double Cave = vC[j] / vNumSamples[j];
        double Csqave = 0.0;
        vC[j] = Cave;
        if (report_rms) {
          assert(L+1 <= vCrms.size());
          Csqave = vCrms[j] / vNumSamples[j];
          vCrms[j] = Csqave - Cave*Cave;
          if (vCrms[j] < 0.0)
            vCrms[j] = 0.0;
          else
            vCrms[j] = sqrt(vCrms[j]);
        }
      }
      else {
        vC[j] = 0.0;
        vCrms[j] = 0.0;
      }
    }
  } //Finalize()

  /// @brief  Sum all of the entries in vC.  Do this after invoking Finalize()
  Scalar
  Integrate() {
    Scalar integral_of_C = 0.0;
    assert(L+1 <= vC.size());
    for (size_t j=0; j <= L; ++j) {
      if ((vNumSamples[j] > 0) && (vC[j] > threshold * vC[0]))
        integral_of_C += vC[j];
      else
        break;
    }
    return integral_of_C;
  }

  /// @brief
  /// Find the j such that vC[j]/vC[0] drops below "thresh".
  /// Use linear interpolation to find the fractional j value close to the
  /// place where the plot of vC[j]/vC[0] drops below that threshold.
  /// If vC[j]/vC[0] remains above the threshold for all j values, return -1.0.
  Scalar
  ThresholdCrossing(Scalar thresh) {
    assert(L+1 <= vC.size());
    size_t j_prev = 0;
    Scalar j_threshold = 0.0;
    Scalar delta_j = -1.0;
    for (size_t j=1; j < vC.size(); j++) {
      if (vNumSamples[j] == 0) //ignore j entries which lack data (if present)
        continue;
      if (vC[j] < thresh * vC[0]) {
        delta_j =
          (thresh*vC[0] - vC[j_prev])/(vC[j]-vC[j_prev]);
        break;
      }
      j_prev = j;
    }
    if (delta_j >= 0.0) {
      Scalar j_thresh = j_prev + delta_j;
      return j_thresh;
    }
    else return -1.0;
  }

  /// @brief  Calculate the correlation length
  ///     (Note: For time series data, this is called the "correlation time".)
  Scalar
  GuessCorrelationLength() {
    // Pick a point along the curve ("j_thresh").
    // Estimate correlation length by observing how much vC[j_thresh] has
    // decayed, and fitting this to an exponential decay.
    Scalar j_thresh = -1.0;
    Scalar C_thresh = -1.0;
    if (persistence_length_threshold > -1.0) {
      // If the "persistence_length_threshold" parameter was specified, then we
      // set "j_thresh" to the point on the curve that crosses this threshold.
      j_thresh = ThresholdCrossing(persistence_length_threshold);
      C_thresh = persistence_length_threshold;
    }

    // --- Optional: If L was specified, set j_thresh = L ---
    //else if (L > 0) {
    //  // If not, then we set "j_thresh" to the last entry in the curve
    //  j_thresh = L;
    //  C_thresh = vC[j_thresh];
    //}
    // In retrospect, this was a bad idea.  I want to allow the user to
    // specify extremely large L values.  In that case, it's probably more
    // accurate to estimate the correlation length using the Integrate() method.


    // Now choose which method to use to calculate the persistence length:
    Scalar persistence_length;

    // If j_thresh and C_thresh are both positive (not pathelogical), then
    // assume the curve is a decaying exponential.  In that case we can
    // use j_thresh and C_thresh to estimate the rate of decay.
    // (The persistence length is one over this rate.)
    if ((j_thresh > 0.0) && (C_thresh > 0.0)) {
      persistence_length = -j_thresh / log(C_thresh / vC[0]);
    }
    else {
      // Otherwise, estimate the correlation length from the
      // integral of the correlation function.
      // (This is numerically unstable, so do this only as a last resort.)
      Scalar integral_of_C = Integrate();
      persistence_length = integral_of_C / vC[0];
    }

    return persistence_length;

  } //NdAutocrr::GuessCorrelationLength()



private:

  void Resize(size_t _L) {
    L = _L;
    vC.resize(L+1);
    vCrms.resize(L+1);
    vNumSamples.resize(L+1);
  }

  /// @brief  Choose the domain of the correlation function
  ///         C(j) is defined from 0 to L-1

  size_t ChooseL(size_t L_single)
  {
    size_t L_backup = L;
    if (is_periodic) {
      if ((L <= 0) || (L > L_single/2))
        L = L_single/2; //default value
    }
    else
    {
      if (L == 0)
        L = L_single/2; // default
      else if (L > L_single-1) 
        L = L_single-1;
    } // if (is_periodic)

    //if ((vC.size() > L+1) && (! is_periodic))
    //  L = vC.size()-1;

    if (L < L_backup)
      L = L_backup;

    if (L+1 > vC.size()) {
      // Allocate enough space to store results from the incomming data.
      size_t size_diff = L+1 - vC.size();
      vC.insert(vC.end(), size_diff, 0.0);
      vNumSamples.insert(vNumSamples.end(), size_diff, 0);
      assert(L+1 == vC.size());
      assert(L+1 == vNumSamples.size());
      if (report_rms) {
        vCrms.insert(vCrms.end(), size_diff, 0.0);
        assert(L+1 == vCrms.size());
      }
    }
    Resize(L);
    return L;
  } //ChooseL

}; //class NdAutocrr





#endif //#ifndef _NDAUTOCRR_HPP
