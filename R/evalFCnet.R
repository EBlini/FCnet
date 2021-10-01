#' Evaluate model's prediction against the observed true values
#'
#' @param model an object returned from a call to cv_FCnet
#' @return A dataframe with model's goodness of fit statistics.
#'
#' @export

evalFCnet= function(model, family){

  # Model performance metrics
  RMSE= sqrt(mean((model$predictions - model$y)^2))

  fit= model$fit$glmnet.fit

  which_lambda= which(fit$lambda== model$lambda)

  tLL= fit$nulldev - fit$nulldev * (1 - fit$dev.ratio)[which_lambda]
  #tLL= deviance(fit)[which_lambda]
  #tLL= sum(dnorm((model$predictions - model$y), log= T))

  k= fit$df[which_lambda]

  n= fit$nobs

  AIC= - 2*tLL + 2 * k

  BIC= log(n) * k - 2*tLL

  #join
  res= data.frame(
    RMSE= RMSE,
    AIC= AIC,
    BIC= BIC
  )

  if(family== "binomial"){

    p= ifelse(model$predictions>0, 1, 0)

    res[["Accuracy"]]= sum(p== model$y)/length(model$y)

  }

  return(res)
}
