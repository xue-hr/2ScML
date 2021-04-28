OneSample_2ScML <- function(Y,D,Z,
                            K1_vec,K2_vec,
                            tau1 = 1e-5,
                            tau2 = 1e-5)
{
  n = length(Y)
  p = ncol(Z)
  
  Y = scale(Y,scale = F)
  D = scale(D,scale = F)
  Z = scale(Z,scale = F)
  
  ### Stage1
  BIC_vec = NULL
  for(K in K1_vec)
  {
    Stage1_Fit = TLC(Y = D,X = Z,K = K,
                     tau = tau1)
    nonzero_ind = which(Stage1_Fit!=0)
    lm_stage1 = lm(D ~ 0 + Z[,nonzero_ind])
    BIC_vec = c(BIC_vec,n*log(sum((lm_stage1$residuals)^2)) + 
                  log(n)*length(nonzero_ind))
  }
  K1 = K1_vec[which.min(BIC_vec)]
  Stage1_Fit = TLC(Y = D,X = Z,K = K1,
                   tau = tau1)
  ind1 = which(Stage1_Fit!=0)
  lm_stage1 = lm(D~0+Z[,ind1])
  Dhat = predict(lm_stage1)
  
  ### Stage2
  BIC_vec = NULL
  for(K in K2_vec)
  {
    Stage2_Fit = TLC(Y = Y,X = cbind(Dhat,Z),K = K,
                     tlc_weight = c(0,rep(1,p)),
                     tau = tau2)
    nonzero_ind = which(Stage2_Fit!=0)
    lm_stage2 = lm(Y ~ 0 + cbind(Dhat,Z)[,nonzero_ind])
    BIC_vec = c(BIC_vec,n*log(sum((lm_stage2$residuals)^2)) + 
                  log(n)*length(nonzero_ind))
  }
  K2 = K2_vec[which.min(BIC_vec)]
  Stage2_Fit = TLC(Y = Y,X = cbind(Dhat,Z),K = K2,
                   tlc_weight = c(0,rep(1,p)),
                   tau = tau2)
  ind2 = which(Stage2_Fit!=0)[-1]-1
  
  ### result
  TLP_Wald_test = OneSample_se(Y,D,Z,
                               ind1,ind2)
  return(list(K1 = K1, K2 = K2, 
              stage1_ind = ind1,
              stage2_ind = ind2,
              beta_est = TLP_Wald_test[1],
              beta_se = TLP_Wald_test[2]))
  
}