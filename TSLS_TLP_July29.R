##################################################
### Two Stage TLP
##################################################
TSLS_TLP <- function(Y,X,sim_bed,K_v_stage1,K_v_stage2,tau)
{
  ### Y is the outcome
  ### X is the exposure
  ### sim_bed are instruments
  ### K_v_stage1 are K's for cv in stage 1
  ### K_v_stage2 are K's for cv in stage 2
  ### tau is the parameter for TLP

  Y = Y
  X = X
  sim_bed = sim_bed
  
  sample_size = length(Y)
  p = ncol(sim_bed)
  
  ### TLP for stage 1
  X_s = scale(X) * sqrt(sample_size / (sample_size-1))
  snp_bed_s = scale(sim_bed) * sqrt(sample_size / (sample_size-1))
  
  rownames(X_s) = NULL
  rownames(snp_bed_s) = NULL
  colnames(snp_bed_s) = NULL
  
  TLP_cv_stage1 = cv_Lasso2_Constrained_TLP(Y = X_s, X = cbind(snp_bed_s),
                                            K_v = K_v_stage1, tau = tau, num_fold = 5,
                                            weight_tlp = c(rep(1,p)),
                                            maxit_tlp = 1000,tol_tlp=1e-7,
                                            maxit_lasso = 2000, tol_lasso = 1e-7,
                                            rho = 1)
  K1 = TLP_cv_stage1$K
  TLP_stage1 = Lasso2_Constrained_TLP(Y = X_s, X = cbind(snp_bed_s),
                                      K = K1, tau = tau,
                                      weight_tlp = c(rep(1,p)),
                                      maxit_tlp = 1000,tol_tlp=1e-7,
                                      maxit_lasso = 2000, tol_lasso = 1e-7,
                                      rho = 1)
  ind = which(TLP_stage1$beta!=0)
  
  
  ### stage1
  lm_stage1 = lm(X~sim_bed[,ind])
  Xhat = predict(lm_stage1)
  
  ### TLP Cross-validation
  Y_s = scale(Y) * sqrt(sample_size / (sample_size-1))
  Xhat_s = scale(Xhat) * sqrt(sample_size / (sample_size-1))
  snp_bed_s = scale(sim_bed) * sqrt(sample_size / (sample_size-1))
  
  rownames(Xhat_s) = NULL
  rownames(snp_bed_s) = NULL
  colnames(snp_bed_s) = NULL
  
  TLP_cv_result_lasso2 = cv_Lasso2_Constrained_TLP(Y = Y_s, X = cbind(Xhat_s,snp_bed_s),
                                                   K_v = K_v_stage2, tau = tau, num_fold = 5,
                                                   weight_tlp = c(0,rep(1,p)),
                                                   maxit_tlp = 1000,tol_tlp=1e-7,
                                                   maxit_lasso = 2000, tol_lasso = 1e-7,
                                                   rho = 1)
  ### Best K and TLP model
  K = TLP_cv_result_lasso2$K
  if(K == 0)
  {
    
    alldata_non_zero = 1
    
    Final_full_model = lm(Y_s ~ cbind(Xhat_s,snp_bed_s)[,alldata_non_zero])
    
    ### LRT
    LRT = 0
    
  } else{
    TLP_all_data = Lasso2_Constrained_TLP(Y = Y_s, X = cbind(Xhat_s,snp_bed_s),
                                          K = K, tau = tau,
                                          weight_tlp = c(0,rep(1,p)),
                                          maxit_tlp = 1000,tol_tlp=1e-7,
                                          maxit_lasso = 2000, tol_lasso = 1e-7,
                                          rho = 1)
    alldata_non_zero = which(TLP_all_data$beta!=0)
    
    Final_full_model = lm(Y_s ~ cbind(Xhat_s,snp_bed_s)[,alldata_non_zero])
    Final_null_model = lm(Y_s ~ snp_bed_s[,alldata_non_zero[-1]-1])
    
    ### LRT
    LRT = sample_size*log(sum((Final_null_model$residuals)^2) / sum((Final_full_model$residuals)^2))
    
  }
  
  return(list(K1 = K1, K2 = K, Stage1_Z = ind, Stage2_Z = alldata_non_zero[-1]-1,
              LRT = LRT, Wald =  summary(Final_full_model)$coefficients[2,1:4]
              ))
}


##################################################
### Calculate Corrected Standard Errors
##################################################
TwoStageStandardError <- function(Y,D,Z,stage1_ind,stage2_ind,n)
{
  lm_stage1 = lm(D~Z[,stage1_ind])
  Xhat = predict(lm_stage1)
  
  lm_True = summary(lm(Y~cbind(Xhat,Z[,stage2_ind])))
  
  X = cbind(Xhat,Z[,stage2_ind])
  inv_Cap_Sigma = solve((t(X)%*%X)/n)
  ZA = Z[,stage1_ind]
  PZA = ZA%*%solve(t(ZA)%*%ZA)%*%t(ZA)
  Cap_Psi = t(X)%*%PZA%*%X/n
  
  Xstar = cbind(D,Z[,stage2_ind])
  hatsigmacombsq = sum((Y - X%*%(lm_True$coefficients[-1,1]))^2)/n
  hatsigma1sq = sum((Y - Xstar%*%(lm_True$coefficients[-1,1]))^2)/n
  v = hatsigmacombsq*inv_Cap_Sigma - 
    (hatsigmacombsq-hatsigma1sq)*((inv_Cap_Sigma%*%Cap_Psi%*%inv_Cap_Sigma))
  return(c(lm_True$coefficients[2,1],sqrt(v[1,1]/n)))
}


