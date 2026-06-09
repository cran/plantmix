// @file GMA_SMA.hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type GMA_SMA(objective_function<Type>* obj)
{
  DATA_VECTOR(y);            // vector of responses
  DATA_MATRIX(X);            // design matrix of fixed effects
  DATA_MATRIX(Z_GMA);        // design matrix of GMAs
  DATA_MATRIX(K);            // correlation matrix for GMAs (kinship)
  DATA_MATRIX(Z_SMA1);       // first design matrix of SMAs, optional
  DATA_MATRIX(Z_SMA2);       // second design matrix of SMAs, optional

  PARAMETER_VECTOR(beta);    // vector of fixed effects
  PARAMETER_VECTOR(GMA);     // vector of random GMAs
  PARAMETER(log_sigma_GMA);  // log std dev of random GMAs
  PARAMETER(log_sigma_e);    // vector of errors
  PARAMETER_VECTOR(SMA1);    // first vector of random SMAs
  PARAMETER(log_sigma_SMA1); // log std dev of random SMAs
  PARAMETER_VECTOR(SMA2);    // second vector of random SMAs
  PARAMETER(log_sigma_SMA2); // log std dev of random SMAs

  Type nll = Type(0.0);

  int n = X.rows();
  vector<Type> m(n);
  m = X * beta;

  m = m + Z_GMA * GMA;
  int q_G = Z_GMA.cols();
  MVNORM_t<Type> mvn_GMA_cor(K);
  vector<Type> vec_sigma_GMA(q_G);
  vec_sigma_GMA.fill(exp(log_sigma_GMA));
  nll += VECSCALE(mvn_GMA_cor, vec_sigma_GMA)(GMA);

  REPORT(beta);
  REPORT(GMA);
  Type var_GMA = exp(2 * log_sigma_GMA);
  REPORT(var_GMA);

  if(Z_SMA1.cols() > 1){
    m = m + Z_SMA1 * SMA1;
    nll -= sum(dnorm(SMA1, 0, exp(log_sigma_SMA1), true));
    REPORT(SMA1);
    Type var_SMA1 = exp(2 * log_sigma_SMA1);
    REPORT(var_SMA1);
  }
  if(Z_SMA2.cols() > 1){
    m = m + Z_SMA2 * SMA2;
    nll -= sum(dnorm(SMA2, 0, exp(log_sigma_SMA2), true));
    REPORT(SMA2);
    Type var_SMA2 = exp(2 * log_sigma_SMA2);
    REPORT(var_SMA2);
  }

  nll -= sum(dnorm(y, m, exp(log_sigma_e), true));
  Type var_err = exp(2 * log_sigma_e);
  REPORT(var_err);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
