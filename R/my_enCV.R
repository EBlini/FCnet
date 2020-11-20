my_enCV= function(y, #dependent variable, typically behavior
                  x, #independent variables, tipycally neural measures
                  alpha= seq(0, 1, by= 0.1),
                  lambda= 10^seq(-5, 5, length.out = 200),
                  nfolds= nrow(x),
                  rep_cv= 1){
  
  #ensure you are working with matrices
  if(class(y)[1]== "data.frame"){y= data.matrix(y)}
  if(class(x)[1]== "data.frame"){x= data.matrix(x)}
  
  
  #to add here: what to do with missing values?
  
  #crossvalidate alpha and/or lambda, but only if vectors are supplied
  if (is.vector(alpha) | is.vector(lambda)){
    
    #crossvalidate n times, take the median of the recovered parameters
    cv_ridge= lapply(1:rep_cv, function(time){
      
      cva.glmnet(x= x, y= y, 
                 alpha = alpha, 
                 lambda = lambda, 
                 type.measure = c("mae"), 
                 nfolds= nfolds,
                 grouped= F,
                 intercept= F,
                 standardize= F
      )
    })
    
    pars= lapply(cv_ridge, function(p)get_model_params(p))
    
    pars= do.call(rbind, pars)
    
    lambda= median(pars$lambdaMin)
    
    alpha= median(pars$alpha)
    
  }
  
  #return best parameters
  return(list(alpha= alpha, lambda= lambda))
  
  
}
