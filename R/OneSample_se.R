OneSample_se <- function(Y,D,Z,stage1_ind,stage2_ind)
{
  n = length(Y)
  lm_stage1 = lm(D~Z[,stage1_ind])
  Xhat = predict(lm_stage1)
  
  lm_stage2 = summary(lm(Y~cbind(Xhat,Z[,stage2_ind])))
  
  X = cbind(Xhat,Z[,stage2_ind])
  inv_Cap_Sigma = solve((t(X)%*%X)/n)
  ZA = Z[,stage1_ind]
  PZA = ZA%*%solve(t(ZA)%*%ZA)%*%t(ZA)
  Cap_Psi = t(X)%*%PZA%*%X/n
  
  Xstar = cbind(D,Z[,stage2_ind])
  hatsigmacombsq = sum((Y - X%*%(lm_stage2$coefficients[-1,1]))^2)/n
  hatsigma1sq = sum((Y - Xstar%*%(lm_stage2$coefficients[-1,1]))^2)/n
  v = hatsigmacombsq*inv_Cap_Sigma - 
    (hatsigmacombsq-hatsigma1sq)*((inv_Cap_Sigma%*%Cap_Psi%*%inv_Cap_Sigma))
  return(c(lm_stage2$coefficients[2,1],sqrt(v[1,1]/n)))
}


