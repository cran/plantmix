/// @file lm.hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type lm(objective_function<Type>* obj)
{
  DATA_VECTOR(y); // n-vector of data
  DATA_MATRIX(X); // n x p design matrix of random effects

  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma);

  Type nll = Type(0.0);

  int n = y.size();
  vector<Type> m(n);
  m = X * beta;
  for(int i=0; i<n; i++){
    nll -= dnorm(y(i), m(i), exp(log_sigma), true);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
