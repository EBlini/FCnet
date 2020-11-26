
#' Internal function for retrieving crossvalidated parameters of cv_FCnet
#'
#' Internal function to retrieve the alpha and lambda parameters
#' resulting in the minimum crossvalidation error.
#'
#' @param fit The default color palette for plotting matrices (a vector of colors).
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

  res= data.frame(
    alpha= alpha[best],
    lambdaMin= lambdaMin[best],
    lambdaSE= lambdaSE[best],
    error= error[best]
  )

  return(res)

}
