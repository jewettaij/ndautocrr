#ifndef CUSTOM_INNER_PRODUCT_H
#define CUSTOM_INNER_PRODUCT_H

#include <vector>
#include <cmath>
using namespace std;

// This function is used to multiply entries together (eg. x(i) and x(i+j))
// whcn computing the correlation function
// C(j) = ⟨(**x**(i)-⟨**x**⟩)⋅(**x**(i+j)-⟨**x**⟩)⟩

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
