#this helper function extracts crossvalidated parameters:
#retrieved: https://stackoverflow.com/questions/54803990/extract-the-best-parameters-from-cva-glmnet-object
get_model_params <- function(fit) {
  alpha <- fit$alpha
  lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
  error <- sapply(fit$modlist, function(mod) {
    min(mod$cvm)
  })
  best <- which.min(error)
  data.frame(
    alpha = alpha[best],
    lambdaMin = lambdaMin[best],
    lambdaSE = lambdaSE[best],
    error = error[best]
  )
}
