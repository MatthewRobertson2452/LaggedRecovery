
// To demonstrate AD. A function that would be difficult to differentiate
// analytically.

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(a);
  PARAMETER(b);
  //PARAMETER(c);
  PARAMETER(log_sigma);

  DATA_VECTOR(global_d);
  DATA_VECTOR(local_d);

  int n = local_d.size();

  //Type a = exp(log_a);
  //Type b = exp(log_b);
  Type sigma = exp(log_sigma);
  vector<Type> pred_d(n);
  vector<Type> resid(n);

  Type nll = 0.0;

  pred_d = a*pow(global_d,b);
  resid = local_d-pred_d;


  nll = -dnorm(resid, 0.0, sigma, true).sum();

  REPORT(a);
  REPORT(b);
  //REPORT(c);
  REPORT(log_sigma);

  return(nll);
}

