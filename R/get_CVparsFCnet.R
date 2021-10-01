#adapted from:
#https://stackoverflow.com/questions/63171921/is-there-a-way-in-r-to-determine-aic-from-cv-glmnet
glmnet_cv_aicc <- function(f= fit$modlist,
                           l= optionsFCnet("whichLambda"),
                           ret= optionsFCnet("cv.criterion")){

  whlm <- sapply(f, function(fa) which(fa$lambda == fa[[l]]))

  sapply(1:length(whlm), function(n){
    with(f[[n]]$glmnet.fit,
         {
           tLL <- nulldev - nulldev * (1 - dev.ratio)[whlm[n]]
           k <- df[whlm[n]]
           n <- nobs
           res= list('AIC' = - 2*tLL + 2 * k,
                     'BIC' = log(n) * k - 2*tLL)
           return(res[[ret]])
         })
  })

}

#' Internal function for retrieving crossvalidated parameters of cv_FCnet
#'
#' Internal function to retrieve the alpha and lambda parameters
#' resulting in the minimum crossvalidation error.
#'
#' @param fit An object returned by `glmnetUtils::cva.glmnet()`.
#'
#' @return The crossvalidated parameters with the associated error

#' @export

get_CVparsFCnet= function(fit){

  alpha= fit$alpha

  lambdaMin= sapply(fit$modlist, `[[`, "lambda.min")

  lambdaSE= sapply(fit$modlist, `[[`, "lambda.1se")

  lambda= ifelse(optionsFCnet("whichLambda")== "lambda.min",
                 lambdaMin,
                 lambdaSE)


  if(optionsFCnet("cv.criterion")== "error"){

    error= sapply(fit$modlist, function(mod){

      as.numeric(mod$cvm[which(mod$lambda== lambda)])

    }) } else {

      error= glmnet_cv_aicc(f= fit$modlist,
                            l= optionsFCnet("whichLambda"),
                            ret= optionsFCnet("cv.criterion"))

    } # end cv.criterion


  best= which.min(error)

  #you can use [[1]] here because Ncomp come in batches
  N_comp= fit$modlist[[1]]$glmnet.fit$dim[1]

  res= data.frame(
    alpha= alpha[best],
    lambda.min= lambdaMin[best],
    lambda.1se= lambdaSE[best],
    error= error[best],
    best= best,
    N_comp= N_comp
  )

  return(res)

}
