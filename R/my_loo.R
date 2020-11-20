my_loo= function(y, #dependent variable, typically behavior
                 x, #independent variables, tipycally neural measures
                 alpha= seq(0, 1, by= 0.1),
                 lambda= 10^seq(-5, 5, length.out = 200),
                 nfolds= nrow(x),
                 rep_cv= 1){
  
  lapply_over= 1:nrow(x)
  
  loo= future_lapply(lapply_over, function(r){
    
    new_x= x[-r,]
    
    new_y= y[-r,]
    
    fit= my_en(y = new_y, x = new_x, 
               alpha = alpha, lambda = lambda,
               nfolds= nfolds, rep_cv= rep_cv)
    
    p= predict(fit$fit, s= fit$lambda, 
               newx= data.matrix(x[r,]))
    
    return(list(prediction= p, alpha= fit$alpha, 
                lambda= fit$lambda, coeffs= fit$coeffs))
    
  }, future.seed= T)
  
  #extract and reshape all relevant information
  prediction= sapply(loo, function(x)x[["prediction"]])
  
  alpha= sapply(loo, function(x)x[["alpha"]])
  
  lambda= sapply(loo, function(x)x[["lambda"]])
  
  coeffs= lapply(loo, function(x){data.frame(x$coeffs)})
  coeffs= do.call(rbind, coeffs)
  coeffs$Subject= rep(lapply_over, 
                      each= nrow(coeffs)/length(lapply_over))
  
  #fit statistics
  pars= eval_results(true = as.vector(unlist(y)), 
                     predicted = prediction)
  
  R2= pars[,"Rsquare"]
  MSE= pars[,"MSE"]
  RMSE= pars[,"RMSE"]
  
  return(list(R2= R2, MSE= MSE, RMSE= RMSE,
              predicted= prediction, alpha= alpha,
              lambda= lambda, coeffs= coeffs))
  
}
