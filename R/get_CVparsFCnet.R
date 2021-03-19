
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

  error= sapply(fit$modlist, function(mod){

    min(mod$cvm)

  })

  best= which.min(error)

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
