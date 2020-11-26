#' Evaluate model's prediction against the observed true values
#'
#' @param true Observed (fitted) values.
#' @param predicted Predicted (model-derived) values.
#' @return A dataframe with model's goodness of fit statistics.
#'
#' @export

evalFCnet= function(true, predicted){

  # Model performance metrics
  RMSE= sqrt(mean((predicted - true)^2))
  MSE= RMSE^2
  Rsquare= as.numeric(cor(true, predicted)^2)

  res= data.frame(
    MSE= MSE,
    RMSE= RMSE,
    R2= Rsquare
  )

  return(res)
}
