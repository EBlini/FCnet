#' Leave-One-Out cross-validated fit of elastic nets for Functional Connectivity data
#'
#' This function is a wrapper around `glmnet::glmnet()` as called from
#' `FCnet::FCnet()`. For extended documentation,
#' the readers are encouraged to consult the original source of `glmnet::glmnet()`
#' and its vignette. `glmnet::glmnet()` fits a robust linear model through
#' penalized maximum likelihood computed via the lasso or elastic net
#' regularization path.
#' `FCnetLOO()` requires two objects at minimum: `y` is a vector or data.frame
#' with exactly one column, corresponding to the (behavioral) score to predict; `x`
#' is a data.frame or a list of lists with an entry named "Weights",
#' which includes the independent variables. `x` can be - and is meant to be -
#' one object created by `reduce_featuresFC()`, but this is not strictly necessary.
#' The best model and hyperparameters are retrieved in inner loops through `cv_FCnet()`;
#' details of the crossvalidation procedures can be passed as arguments to
#' `FCnetLOO()` if necessary.
#' A call to `FCnetLOO()` returns a list including goodness of fit measures
#' for the outer loop, a data.frame including coefficients for all nested models,
#' vectors of crossvalidated parameters, and predicted scores.
#' Note that dependent variables are scaled by default through `scale_y`.
#' The `ParallelLOO` option is recommended for speed. In order to use parallel
#' computing, `future.apply` must be installed, your machine should have
#' multiple cores available, and parallel computing should be prompted by
#' the user (e.g. via `plan(multisession`).
#'

#'
#'
#' @param y The dependent variable, typically behavioral scores to predict.
#' This can be a vector or a single data.frame column.
#' @param x The independent variables, typically neural measures that have
#' been already summarised through data reduction techniques
#' (e.g. ICA, PCA): an object created by `reduce_featuresFC()` will do. If such
#' an object is passed to this function, the "Weights" slot is taken as x.
#' A list can be passed to this function: in this case the function needs an
#' entry named "Weights". Otherwise, a data.frame can be passed to x.
#' @param alpha Value(s) that bias the elastic net toward ridge regression
#' (alpha== 0) or LASSO regression (alpha== 1). If a vector of alpha values
#' is supplied, the value is optimized through crossvalidation.
#' It defaults to a vector ranging from 0 to 1 with steps of 0.1.
#' The crossvalidated alpha is returned.
#' @param lambda Regularization parameter for the regression,
#' see `glmnet::glmnet()`. Lambda must be a vector with length>1.
#' When a vector of lambda values is supplied, the value of lambda
#' is optimized through internal crossvalidation. It defaults to a vector
#' ranging from 10^-5 to 10^5 with 200 values in logarithmic steps.
#' The crossvalidated optimal lambda is returned.
#' @param parallelLOO If TRUE - recommended, but not the default - uses
#' `future.apply::future_lapply()` for the outer loops: `future.apply` must be
#' installed, the machine should have multiple cores available for use,
#' and threads should be defined explicitly beforehand by the user
#' (e.g. by calling `plan(multisession)`).
#' @param scale_y Whether y should be scaled prior to fit. Default, TRUE, scales y with
#' `scale()`.
#' @param cv.type.measure The measure to minimize in crossvalidation inner loops.
#' Differently from `glmnetUtils::cva.glmnet()` the deafult is the mean absolute error.
#' @param intercept whether to fit (TRUE) or not (FALSE) an intercept to the model.
#' @param standardize Whether x must be standardized.
#' `glmnet::glmnet()` the default is FALSE as we assume predictors are already either
#' summarised with PCA or ICA (and therefore scaled) or drawn from normalized FC matrices.
#' @param ... Other parameters passed to `glmnetUtils::cva.glmnet()` or
#' `glmnet::glmnet()`.
#'
#' @return A model produced by `glmnet::glmnet()`.
#' Crossvalidated best alpha and lambda values.
#' If requested: a vector of observed (y) and predicted (predicted) values;
#' model's coefficients; goodness of fit statistics.
#'
#' @export


FCnetLOO= function(y,
                   x,
                   alpha= seq(0, 1, by= 0.1),
                   lambda= rev(10^seq(-5, 5, length.out = 200)),
                   parallelLOO= F,
                   scale_y= T,
                   type.measure= optionsFCnet("cv.type.measure"),
                   intercept= optionsFCnet("intercept"),
                   standardize= optionsFCnet("standardize"),
                   ...){

  #ensure you are working with matrices
  y= data.matrix(y)

  #scaling if requested
  if(scale_y){y= scale(y)}

  #further reformatting of x
  if(class(x)[1]== "list"){

    x= data.matrix(x$Weights)

  } else {

    if(class(x)[1]== "data.frame"){x= data.matrix(x)}

  }


  #indices for LOO
  lapply_over= 1:nrow(x)

  #main sequence here: parallel or not

  if(parallelLOO== TRUE){

    loo= future.apply::future_lapply(lapply_over, function(r){

      new_x= x[-r,]

      new_y= y[-r,]

      fit= FCnet(y = new_y,
                 x = new_x,
                 alpha = alpha,
                 lambda = lambda,
                 fit_stats = T,
                 type.measure= type.measure,
                 intercept= intercept,
                 standardize= standardize,
                 ...
                 )

      p= predict(fit$fit,
                 s= fit$lambda,
                 newx= t(data.matrix((x[r,]))))

      p= as.numeric(p)

      return(list(prediction= p,
                  alpha= fit$alpha,
                  lambda= fit$lambda,
                  coeffs= fit$coeffs))

    }, future.seed = T)

  } else {

    loo= lapply(lapply_over, function(r){

      new_x= x[-r,]

      new_y= y[-r,]

      fit= FCnet(y = new_y,
                 x = new_x,
                 alpha = alpha,
                 lambda = lambda,
                 fit_stats = T,
                 type.measure= type.measure,
                 intercept= intercept,
                 standardize= standardize,
                 ...)

      p= predict(fit$fit,
                 s= fit$lambda,
                 newx= t(data.matrix((x[r,])))
                 )

      p= as.numeric(p)

      return(list(prediction= p,
                  alpha= fit$alpha,
                  lambda= fit$lambda,
                  coeffs= fit$coeffs))

    })

  } #end if parallelLOO

  #extract and reshape all relevant information
  prediction= sapply(loo, function(x)x[["prediction"]])

  alpha= sapply(loo, function(x)x[["alpha"]])

  lambda= sapply(loo, function(x)x[["lambda"]])

  coeffs= lapply(loo, function(x){data.frame(x$coeffs)})
  coeffs= do.call(rbind, coeffs)
  coeffs$ID= rep(lapply_over,
                 each= nrow(coeffs)/length(lapply_over))

  #fit statistics
  pars= evalFCnet(true = as.vector(unlist(y)),
                  predicted = prediction)

  R2= pars[,"R2"]
  MSE= pars[,"MSE"]
  RMSE= pars[,"RMSE"]

  #wrap-up info
  res= list(R2= R2,
            MSE= MSE,
            RMSE= RMSE,
            predicted= prediction,
            alpha= alpha,
            lambda= lambda,
            coeffs= coeffs,
            y= y)

  return(res)

}
