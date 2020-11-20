eval_results= function(true, predicted){
  # Model performance metrics
  data.frame(
    MSE = RMSE(obs = true, pred = predicted)^2,
    RMSE = RMSE(obs = true, pred = predicted),
    Rsquare = R2(obs = true, pred = predicted)
  )
}