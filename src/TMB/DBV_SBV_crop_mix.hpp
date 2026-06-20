/// @file DBV_SBV_crop_mix.hpp

// Model:
// * intercrops: $Y_IC = X_IC B_IC + Z_{DS} BV_f + E_IC$
// * sole crops:
//     * focal species: $y_{SC_f} = X_{SC_f} \beta_{SC_f} + Z_D DBV_f + Z_D SIGV_f + e_{SC_f}$
//     * tester species: $y_{SC_t} = X_{SC_t} \beta_{SC_t} + e_{SC_t}$

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type DBV_SBV_crop_mix(objective_function<Type>* obj)
{
  // Caution, TMB does not allows DATA or PARAMETER in if statements!

  // IC
  DATA_ARRAY(Y_IC);
  DATA_SPARSE_MATRIX(X_IC);
  DATA_SPARSE_MATRIX(Z_DS_f);
  DATA_SPARSE_MATRIX(Z_DxS_f);
  DATA_MATRIX(K);
  DATA_MATRIX(Kmixpair);
  DATA_SPARSE_MATRIX(invKmixpair);
  // SC focal
  DATA_VECTOR(y_SC_f);
  DATA_SPARSE_MATRIX(X_SC_f);
  DATA_SPARSE_MATRIX(Z_D_f);
  // SC tester
  DATA_VECTOR(y_SC_t);
  DATA_SPARSE_MATRIX(X_SC_t);

  // IC
  PARAMETER_MATRIX(B_IC);
  PARAMETER_ARRAY(BV_f);
  PARAMETER_VECTOR(log_sd_BV_f);
  PARAMETER_VECTOR(unconstr_cor_DS_f);
  PARAMETER_ARRAY(DBVxSBV);
  PARAMETER_VECTOR(log_sd_DxS);
  PARAMETER_VECTOR(unconstr_cor_DxS);
  PARAMETER_VECTOR(log_sd_E_IC);
  PARAMETER_VECTOR(unconstr_cor_E_IC);
  // SC focal
  PARAMETER_VECTOR(beta_SC_f);
  PARAMETER_VECTOR(SIGV_f);
  PARAMETER(log_sd_SIGV_f);
  PARAMETER(log_sd_e_SC_f);
  // SC tester
  PARAMETER_VECTOR(beta_SC_t);
  PARAMETER(log_sd_e_SC_t);

  // Importantly, the DBV_f used to model SC data are the same as the ones
  // used to model IC data; they initially are in the input BV_f matrix.
  vector<Type> DBV_f = BV_f.matrix().col(0);

  Type nll = Type(0.0); // negative log-likelihood

  // contrib of BV_f from the focal species to the nll
  UNSTRUCTURED_CORR_t<Type> mvn_u_bv(unconstr_cor_DS_f);
  vector<Type> sd_BV_f = exp(log_sd_BV_f);
  VECSCALE_t<UNSTRUCTURED_CORR_t<Type>> f_bv = VECSCALE(mvn_u_bv, sd_BV_f);
  MVNORM_t<Type> g_bv(K);
  SEPARABLE_t<VECSCALE_t<UNSTRUCTURED_CORR_t<Type>>, MVNORM_t<Type>> h_BV(f_bv, g_bv);
  nll += h_BV(BV_f);
  // simulate BV_f
  array<Type> BV_f_sim(BV_f.rows(), 2);
  SIMULATE {
    h_BV.simulate(BV_f_sim);
    // REPORT(BV_f_sim);
  }

  if(y_SC_f.size() > 1){
    // contrib of SIGV from the focal species to the nll
    MVNORM_t<Type> mvn_SIGV_cor(K);
    vector<Type> sd_SIGV_f(K.cols());
    sd_SIGV_f.fill(exp(log_sd_SIGV_f));
    nll += VECSCALE(mvn_SIGV_cor, sd_SIGV_f)(SIGV_f);
    // simulate SIGV_f
    vector<Type> SIGV_f_sim(SIGV_f.size());
    SIMULATE {
      VECSCALE_t<MVNORM_t<Type>> h_SIGV = VECSCALE(mvn_SIGV_cor, sd_SIGV_f);
      h_SIGV.simulate(SIGV_f_sim);
      // REPORT(SIGV_f_sim);
    }

    // contrib of the SC observations from the focal species to the nll
    vector<Type> m_SC_f(y_SC_f.size());
    m_SC_f = X_SC_f * beta_SC_f + Z_D_f * DBV_f + Z_D_f * SIGV_f;
    SparseMatrix<Type> Id_SC_f(y_SC_f.size(), y_SC_f.size());
    fillSparseMatId(Id_SC_f);
    MVNORM_t<Type> mvn_y_SC_f_cor(Id_SC_f);
    vector<Type> sd_e_SC_f(y_SC_f.size());
    sd_e_SC_f.fill(exp(log_sd_e_SC_f));
    nll += VECSCALE(mvn_y_SC_f_cor, sd_e_SC_f)(y_SC_f - m_SC_f);
    // simulate
    vector<Type> e_SC_f_sim(y_SC_f.size());
    vector<Type> DBV_f_sim = BV_f_sim.matrix().col(0);
    vector<Type> y_SC_f_sim(y_SC_f.size());
    SIMULATE {
      VECSCALE_t<MVNORM_t<Type>> h_e_SC_f = VECSCALE(mvn_y_SC_f_cor, sd_e_SC_f);
      h_e_SC_f.simulate(e_SC_f_sim);
      y_SC_f_sim  =  X_SC_f * beta_SC_f + Z_D_f * DBV_f_sim +
        Z_D_f * SIGV_f_sim + e_SC_f_sim;
      REPORT(y_SC_f_sim);
    }
  }

  if(y_SC_t.size() > 1){
    // contrib of the SC observations from the tester species to the nll
    vector<Type> m_SC_t(y_SC_t.size());
    m_SC_t = X_SC_t * beta_SC_t;
    SparseMatrix<Type> Id_SC_t(y_SC_t.size(), y_SC_t.size());
    fillSparseMatId(Id_SC_t);
    MVNORM_t<Type> mvn_y_SC_t_cor(Id_SC_t);
    vector<Type> sd_e_SC_t(y_SC_t.size());
    sd_e_SC_t.fill(exp(log_sd_e_SC_t));
    nll += VECSCALE(mvn_y_SC_t_cor, sd_e_SC_t)(y_SC_t - m_SC_t);
    // simulate
    vector<Type> e_SC_t_sim(y_SC_t.size());
    vector<Type> y_SC_t_sim(y_SC_t.size());
    SIMULATE {
      VECSCALE_t<MVNORM_t<Type>> h_e_SC_t = VECSCALE(mvn_y_SC_t_cor, sd_e_SC_t);
      h_e_SC_t.simulate(e_SC_t_sim);
      y_SC_t_sim  = X_SC_t * beta_SC_t + e_SC_t_sim;
      REPORT(y_SC_t_sim); 
    }
  }

  // contrib of DBVxSBV to the nll
  UNSTRUCTURED_CORR_t<Type> mvn_u_DxS(unconstr_cor_DxS);
  array<Type> DBVxSBV_sim(DBVxSBV.rows(), DBVxSBV.cols()); // must be in global scope to be accessible
  if(Z_DxS_f.cols() > 2){
    vector<Type> sd_DxS = exp(log_sd_DxS);
    VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > f_DxS = VECSCALE(mvn_u_DxS, sd_DxS);
    if(Kmixpair.cols() > 2){
      MVNORM_t<Type> g_DxS(Kmixpair);
      SEPARABLE_t< VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > , MVNORM_t<Type> > h_DxS(f_DxS, g_DxS);
      nll += h_DxS(DBVxSBV);
      SIMULATE {
        h_DxS.simulate(DBVxSBV_sim);
        // REPORT(DBVxSBV_sim);
      }
    } else if (invKmixpair.cols() > 2){
      GMRF_t<Type> g_DxS(invKmixpair);
      SEPARABLE_t< VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > , GMRF_t<Type> > h_DxS(f_DxS, g_DxS);
      nll += h_DxS(DBVxSBV); 
      SIMULATE {
        h_DxS.simulate(DBVxSBV_sim);
        // REPORT(DBVxSBV_sim);
      }
    }
  }

  // contrib of the IC observations to the nll
  UNSTRUCTURED_CORR_t<Type> mvn_Y(unconstr_cor_E_IC);
  vector<Type> sd_E = exp(log_sd_E_IC);
  VECSCALE_t<UNSTRUCTURED_CORR_t<Type>> f_Y = VECSCALE(mvn_Y, sd_E);
  SparseMatrix<Type> Id_n(Y_IC.rows(), Y_IC.rows());
  fillSparseMatId(Id_n);
  GMRF_t<Type> g_Y(Id_n);
  SEPARABLE_t<VECSCALE_t<UNSTRUCTURED_CORR_t<Type>>, GMRF_t<Type>> h_Y(f_Y, g_Y);
  matrix<Type> M(Y_IC.rows(), Y_IC.cols());
  M = X_IC * B_IC + Z_DS_f * BV_f.matrix();
  if(Z_DxS_f.cols() > 2){
    M = M + Z_DxS_f * DBVxSBV.matrix();
  }
  nll += h_Y(Y_IC - M.vec());
  // simulate
  matrix<Type> M_sim(Y_IC.rows(), Y_IC.cols());
  array<Type> E_IC_sim(Y_IC.rows(), Y_IC.cols());
  array<Type> Y_IC_sim(Y_IC.rows(), Y_IC.cols());
  SIMULATE {
    M_sim = X_IC * B_IC + Z_DS_f * BV_f_sim.matrix();
    if(Z_DxS_f.cols() > 2){
      M_sim = M_sim + Z_DxS_f * DBVxSBV_sim.matrix();
    }
    h_Y.simulate(E_IC_sim);
    for (int i = 0; i < Y_IC_sim.rows(); i++) {
      for (int j = 0; j < Y_IC_sim.cols(); j++) {
        Y_IC_sim(i, j) = M_sim(i, j) + E_IC_sim.matrix()(i, j);
      }
    }
    REPORT(Y_IC_sim);
  }

  // reports ("AD" to get std.err)
  ADREPORT(BV_f);
  vector<Type> BV_IC_f(BV_f.rows());
  BV_IC_f = BV_f.col(0) + BV_f.col(1);
  ADREPORT(BV_IC_f);

  vector<Type> vars_BV_f = exp(2 * log_sd_BV_f);
  ADREPORT(vars_BV_f);

  Type cor_BV_f = Type(0.0);
  cor_BV_f = mvn_u_bv.cov()(0,1);
  ADREPORT(cor_BV_f);

  if(Z_DxS_f.cols() > 2){
    matrix<Type> Cor_DxS(2, 2);
    Cor_DxS = mvn_u_DxS.cov();
    ADREPORT(Cor_DxS);
  }

  vector<Type> vars_E_IC = exp(2 * log_sd_E_IC);
  ADREPORT(vars_E_IC);
  
  matrix<Type> Cor_E_IC(Y_IC.cols(), Y_IC.cols());
  Cor_E_IC = mvn_Y.cov();
  ADREPORT(Cor_E_IC);
  
  ADREPORT(B_IC);
  if(y_SC_f.size() > 1) {
    ADREPORT(beta_SC_f);
    Type var_SIGV_f = exp(2 * log_sd_SIGV_f);
    ADREPORT(var_SIGV_f);
    ADREPORT(SIGV_f);
    Type var_err_SC_f = exp(2 * log_sd_e_SC_f);
    ADREPORT(var_err_SC_f);
    vector<Type> BV_SC_f(BV_f.rows());
    for (int i = 0; i < BV_SC_f.size(); i++)
      BV_SC_f(i) = BV_f(i,0) + SIGV_f(i);
    ADREPORT(BV_SC_f);
  }
  if(y_SC_f.size() > 1) {
    ADREPORT(beta_SC_t);
    Type var_err_SC_t = exp(2 * log_sd_e_SC_t);
    ADREPORT(var_err_SC_t);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
