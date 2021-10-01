# Variable, global to package's namespace.
# This function is not exported to user space and does not need to be documented.
opt= settings::options_manager(
  family= "gaussian",
  cv.type.measure= "mse",
  cv.criterion= "error",
  intercept= F,
  standardize= F,
  thresh= 10^-4,
  whichLambda= "lambda.min",
  consensus_function= median,
  nested= T,

  colorPaletteDefault= c("darkslateblue", "royalblue4",
                       "royalblue1", "cyan2",
                       "white", "white",
                       "darkgoldenrod1", "darkgoldenrod3",
                       "red3", "firebrick4"),
  colorNuances= 20

  )

# User function that gets exported:

#' Set or get options for FCnet
#'
#' @param family Defaults to "gaussian." Experimental support for "binomial" on the way.
#' @param cv.type.measure The measure to minimize in crossvalidation inner loops.
#' The deafult is the mean squared error.
#' @param cv.criterion Whether crossvalidated models are chosen on the basis of
#' their error or - relevant when different number of components are compared -
#' AIC or BIC criteria.
#' @param intercept whether to fit (TRUE) or not (FALSE) an intercept to the model.
#' Useful if y is not scaled and centered.
#' @param standardize Whether x must be standardized internally to `glmnet::glmnet()`.
#' @param whichLambda During crossvalidation, either selects "lambda.min" - which
#' chooses the value of lambda which best minimizes the mean crossvalidation error - or
#' "lambda.1se" - which "which gives the most regularized model such that error is
#' within one standard error of the minimum". The latter would be preferable
#'  for stability.
#' @param thresh Threshold for glmnet to stop converging to the solution.
#'  This parameter can be adapted (e.g. for speed / accuracy tradeoff).
#' @param consensus_function which function is used to create consensus between
#'  coefficients from different models.
#' @param nested Whether crossvalidation is to be performed in a nested fashion.
#' @param colorPaletteDefault The default color palette for plotting matrices (a vector of colors).
#' @param colorNuances Number of nuances along the provided colorPalette.

#' @export

optionsFCnet <- function(...){

  # protect against the use of reserved words.
  settings::stop_if_reserved(...)
  # return/change
  opt(...)

}
