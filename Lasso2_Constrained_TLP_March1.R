library(lasso2)
##################################################
### Lasso2 Constrained TLP
##################################################
Lasso2_Constrained_TLP <- function(Y,X,K,tau,weight_tlp = rep(1,ncol(X)),
                                   maxit_tlp = 1000,tol_tlp=1e-7,
                                   maxit_lasso = 2000, tol_lasso = 1e-7,
                                   rho = 1)
{
  
  ### For the TLP problem, using DC algorithm to iteratively relax
  ### non-convex optimization problem to convex problems
  ### Y is the length n response vector, scaled
  ### X is the n by p design matrix, scaled
  ### K is the constraint parameter
  ### tau is tau in TLP
  ### weight_tlp are weights for variables, if we do not want to
  ### constrain Xi then weight_tlp[i] = 0, otherwise is 1
  ### maxit's are the maximum number of iterations for TLP and lasso respectively
  ### same for tol's
  ### rho is the penalty parameter in ADMM
  
  n = length(Y)
  p = ncol(X)
  
  # Using constrained lasso as initial estimate of beta
  weight_lasso = weight_tlp / tau
  #weight_lasso = weight_tlp
  
  t = K
  ###
  X1 = X
  weight1 = weight_lasso
  zero_ind = which(weight1==0)
  weight1[zero_ind] = 1
  
  X1 = t(t(X1) / weight1)
  
  data1 = as.data.frame(cbind(Y,X1))
  fit_formula = as.formula(paste(colnames(data1)[1], "~."))
  if(length(zero_ind > 0))
  {
    sweep_formula = as.formula(paste("~1+",
                                     paste(colnames(data1)[zero_ind+1],
                                           collapse = "+"),sep=""))
    
  } else {
    sweep_formula = as.formula(paste("~1"))
  }
  
  lasso2fit = l1ce(fit_formula,data1,bound=t,absolute.t = T,
                   sweep.out = sweep_formula,standardize = F)
  beta = lasso2fit$coefficients[-1]
  beta = beta / weight1
  names(beta) = NULL
  ###

  loss_old = sum((Y - X%*%beta)^2)
  loss_new = loss_old
  
  ite_ind = 1
  while(ite_ind <= maxit_tlp)
  {
    #cat(ite_ind,"\n")
    ite_ind = ite_ind+1
    weight_lasso = (weight_tlp / tau) * (abs(beta) <= tau)
    t = K - sum(weight_tlp * (abs(beta) > tau))
    #cat("t = ",t,"\n")
    if(t<=0)
    {
      weight1 = weight_lasso
      zero_ind = which(weight1==0)
      lm1 = lm(Y ~ X[,zero_ind])
      beta1 = (lm1$coefficients)[-1]
      names(beta1) = NULL
      beta = rep(0,p)
      beta[zero_ind] = beta1
      #break()
    } else {
      ###
      X1 = X
      weight1 = weight_lasso
      zero_ind = which(weight1==0)
      weight1[zero_ind] = 1
      
      X1 = t(t(X1) / weight1)
      
      data1 = as.data.frame(cbind(Y,X1))
      if(length(zero_ind > 0))
      {
        sweep_formula = as.formula(paste("~1+",
                                         paste(colnames(data1)[zero_ind+1],
                                               collapse = "+"),sep=""))
        
      } else {
        sweep_formula = as.formula(paste("~1"))
      }
      lasso2fit = l1ce(fit_formula,data1,bound=t,absolute.t = T,
                       sweep.out = sweep_formula,standardize = F)
      beta = lasso2fit$coefficients[-1]
      beta = beta / weight1
      names(beta) = NULL
      ###
    }
    
    #cat(loss_new,"\n")
    #cat(beta,"\n")
    
    loss_new = sum((Y - X%*%beta)^2)
    if(abs(loss_new - loss_old) < tol_tlp)
    {
      break()
    } else{
      loss_old = loss_new
    }
  }
  return(list(beta = beta,loss = loss_new))
}


##################################################
### Lasso2 Constrained TLP Cross Validation
##################################################
cv_Lasso2_Constrained_TLP <- function(Y,X,K_v,tau,num_fold = 5,
                                      weight_tlp = rep(1,ncol(X)),
                                      maxit_tlp = 1000,tol_tlp=1e-7,
                                      maxit_lasso = 2000, tol_lasso = 1e-7,
                                      rho = 1)
{
  n = nrow(X)
  sep_fold = seq(1,n,ceiling(n/num_fold))
  fold_start_end = cbind(sep_fold,c(sep_fold[-1]-1,n))
  
  cv_error_vec = rep(0,length(K_v))
  
  for(i in 1:length(K_v))
  {
    #cat(i,"\n") ### change on March 1
    K = K_v[i]
    cv_error_i = 0
    for(fold_ind in 1:num_fold)
    {
      #cat("Fold = ",fold_ind,"\n") ### change on March 1
      start_ind = fold_start_end[fold_ind,1]
      end_ind = fold_start_end[fold_ind,2]
      test_range = (start_ind:end_ind)
      
      Y_test = Y[test_range]
      X_test = X[test_range,]
      n_test = nrow(X_test)
      Y_test = scale(Y_test) * sqrt(n_test / (n_test-1))
      X_test = scale(X_test) * sqrt(n_test / (n_test-1))
      
      Y_train = Y[-test_range]
      X_train = X[-test_range,]
      n_train = nrow(X_train)
      Y_train = scale(Y_train) * sqrt(n_train / (n_train - 1))
      X_train = scale(X_train) * sqrt(n_train / (n_train - 1))
      if(K == 0)
      {
        
        non_zero = 1
        lm1 = lm(Y_train ~ X_train[,non_zero])
        cv_error_i = cv_error_i + sum((Y_test - cbind(1,X_test[,non_zero])
                                       %*%(lm1$coefficients))^2)
        
      } else{
        TLP_full = Lasso2_Constrained_TLP(Y = Y_train,X = X_train,
                                          K,tau,
                                          weight_tlp = weight_tlp,
                                          maxit_tlp = maxit_tlp, tol_tlp = tol_tlp,
                                          maxit_lasso = maxit_lasso, tol_lasso = tol_lasso,
                                          rho = 1)
        non_zero = which(TLP_full$beta!=0)
        lm1 = lm(Y_train ~ X_train[,non_zero])
        cv_error_i = cv_error_i + sum((Y_test - cbind(1,X_test[,non_zero])
                                       %*%(lm1$coefficients))^2)
        
      }
    }
    cv_error_vec[i] = cv_error_i
  }
  return(list(cv_error_vec = cv_error_vec , K = K_v[which.min(cv_error_vec)]))
}
  