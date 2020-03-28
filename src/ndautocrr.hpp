///   @file ndautocrr.hpp
///   @brief  Calculate the auto-correlation function from multiple time series.
///   @date 2007-4-13

#ifndef _NDAUTOCRR_HPP
#define _NDAUTOCRR_HPP

#include <vector>
#include <cassert>
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
    if (L > 0) {
      // If the user manually specified a domainwidth using the \"-L\" argument
      // then do not use thresholding to make the computation faster
      // (but do use it to calculate the persistence length).
      persistence_length_threshold = threshold;
      threshold = -1.01;
      vC.resize(L+1);
      if (report_rms)
        vCrms.resize(L+1);
      vNumSamples.resize(L+1);
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

      #pragma omp parallel
      {
        size_t jmax = _vvX_id.size();
        if (jmax > L)
          jmax = L;
        #pragma omp for collapse(1)
        for (size_t j=0; j <= jmax; ++j)
        {
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
            L = j-1;  // This will get us out of the loop
          }
        } //for (size_t j=0; j <= L; ++j)
      } //#pragma omp parallel

    } //else clause for "if (is_periodic)"
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
    for (size_t j=0; j < L; ++j) {
      if ((vNumSamples[j] > 0) && (vC[j] > persistence_length_threshold * vC[0]))
        integral_of_C += vC[j];
      else
        break;
    }
    return integral_of_C;
  }

  /// @brief  Calculate the correlation length
  ///     (Note: For time series data, this is called the "correlation time".)
  Scalar
  CorrelationLength() {
    Scalar integral_of_C = Integrate();
    return integral_of_C / vC[0];
  }


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
