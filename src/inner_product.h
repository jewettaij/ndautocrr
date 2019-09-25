#ifndef CUSTOM_INNER_PRODUCT_H
#define CUSTOM_INNER_PRODUCT_H

#include <vector>
#include <cassert>
#include <cmath>
using namespace std;

// Feel free to customize the definition of the inner_product() below.

inline double 
inner_product(const vector<double> &vXa_d,
	      const vector<double> &vXb_d)
{
  long D = vXa_d.size();
  double total = 0.0;
  for (long d=0; d < D; ++d)
    total += vXa_d[d] * vXb_d[d];
  return total;
}


#endif //#ifndef CUSTOM_INNER_PRODUCT_H
