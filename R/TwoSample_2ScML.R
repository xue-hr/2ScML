TwoSample_2ScML <- function(Y,Z,stage1_ind,gamma_hat,
                            Theta0_hat,sigma2_hat,n2,
                            K2_vec,tau2 = 1e-5)
{
  
  n = length(Y)
  p = ncol(Z)
  
  Y = scale(Y,scale = F)
  Z = scale(Z,scale = F)
  
  Dhat = Z[,stage1_ind]%*%gamma_hat
  
  ### Stage2
  BIC_vec = NULL
  for(K in K2_vec)
  {
    Stage2_Fit = TLC(Y = Y,X = cbind(Dhat,Z),K = K,
                     tlc_weight = c(0,rep(1,p)),
                     tau = tau2)
    nonzero_ind = which(Stage2_Fit!=0)
    lm_stage2 = lm(Y ~ 0 + cbind(Dhat,Z)[,nonzero_ind])
    BIC_vec = c(BIC_vec,
                n*log(sum((lm_stage2$residuals)^2)) + log(n)*length(nonzero_ind))
  }
  K2 = K2_vec[which.min(BIC_vec)]
  Stage2_Fit = TLC(Y = Y,X = cbind(Dhat,Z),K = K2,
                   tlc_weight = c(0,rep(1,p)),
                   tau = tau2)
  ind2 = which(Stage2_Fit!=0)[-1]-1
  
  ### result
  TLP_Wald_test = 
    TwoSampleEstSE(Y,Z,stage1_ind,gamma_hat,
                   Theta0_hat,sigma2_hat,n2,
                   ind2)
  
  return(list(K2 = K2, 
              stage2_ind = ind2,
              beta_est = TLP_Wald_test[1],
              beta_se = TLP_Wald_test[2]))
  
}