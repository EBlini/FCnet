my_en= function(y, #dependent variable, typically behavior
                x, #independent variables, tipycally neural measures
                alpha= seq(0, 1, by= 0.1),
                lambda= 10^seq(-5, 5, length.out = 200),
                nfolds= nrow(x),
                rep_cv= 1){
  
  
  #ensure you are working with matrices
  if(class(y)[1]== "data.frame"){y= data.matrix(y)}
  if(class(x)[1]== "data.frame"){x= data.matrix(x)}
  
  pars= my_enCV(y= y, x= x, 
                alpha= alpha, lambda= lambda,
                nfolds= nfolds, rep_cv = rep_cv)
  
  alpha= pars[["alpha"]]
  
  lambda= pars[["lambda"]]
  
  #fit with the final parameters
  fit= glmnet(x,y, alpha= alpha, lambda= lambda, 
              family= 'gaussian',
              standardize= F,
              intercept= F)
  
  #all vars
  all_rn= rownames(coef(fit))
  #rownames with values
  rnwv= coef(fit)[,1] %in% coef(fit)@x
  
  cf= rep(0, length= length(all_rn))
  
  cf[rnwv]= coef(fit)@x
  
  #so by now we have all coefficients 
  #regardless of if they are dropped (e.g. LASSO) 
  coeffs= data.frame(Feature= seq(1, length(all_rn)-1, 1),
                     Coefficient= cf[-1])
  
  p= predict(fit, s= lambda, 
             newx= as.matrix(x))
  
  return(list(fit= fit, 
              alpha= alpha, lambda= lambda, 
              y= y, predicted= p,
              coeffs= coeffs))
  
}
