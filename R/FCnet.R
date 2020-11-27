#' Fit an elastic net model with crossvalidated hyperparameters
#'
#' This function is a wrapper around `glmnet::glmnet()`. For extended documentation,
#' the readers are encouraged to consult the original source.
#' `glmnet::glmnet()` fits a robust linear model through penalized maximum likelihood
#' computed via the lasso or elastic net regularization path.
#' `FCnet()` requires two objects at minimum: `y` is a vector or data.frame
#' with exactly one column, corresponding to the (behavioral) score to predict; `x`
#' is a data.frame or a list of lists with an entry named "Weights",
#' which includes the independent variables. `x` can be - and is meant to be -
#' one object created by `reduce_featuresFC()`, but this is not strictly necessary.
#' The final model finds the best hyperparameters through `cv_FCnet()`; details of
#' the crossvalidation procedures can be passed as arguments to `FCnet()` if
#' necessary. Use of this function is not encouraged because crossvalidation and
#' predictions are made on overlapping sets of observations. A better approach
#' - though computationally more demanding - would be to perform **nested** inner
#' and outer loop, e.g.- via the function `FCnetLOO()`.
#' This is the reason why, in order to evaluate the results of a `FCnet()` call
#' in terms of goodness of fit, such statistics must be required explicitly by
#' setting `fit_stats= T` (sure such statistics can be calculated post-hoc
#' anyway).
#' A call to `FCnet()` returns a `glmnet::glmnet()` model which can be explored
#' further if necessary.
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
#' @param nfolds Number of folds to be created in the crossvalidation of alpha
#' and/or lambda. It defaults to `nrow(x)`, that is Leave-One-Out crossvalidation.
#' It can be set to whatever integer >3 and < nrow(x).
#' @param rep_cv Number of times the crossvalidation procedure must be repeated.
#' It defaults to 1 (as the LOO procedure is deterministic, thus values larger
#' than 1 are redundant). If rep_cv>1, the crossvalidation is repeated rep_cv
#' times, and the median crossvalidated alpha and lambda values across all rep_cv
#' are returned. Useful in order to decrease randomness when kfold crossvalidation
#' is required. Setting this parameter too high though will result in a mere
#' approximation of the LOO.
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

FCnet= function(y,
                x,
                alpha= seq(0, 1, by= 0.1),
                lambda= rev(10^seq(-5, 5, length.out = 200)),
                nfolds= nrow(x),
                rep_cv= 1,
                fit_stats= F,
                type.measure= optionsFCnet("cv.type.measure"),
                intercept= optionsFCnet("intercept"),
                standardize= optionsFCnet("standardize"),
                ...){


  #ensure you are working with matrices
  y= data.matrix(y)

  if(class(x)[1]== "list"){

    x= data.matrix(x$Weights)

  } else {

    if(class(x)[1]== "data.frame"){x= data.matrix(x)}

  }

  pars= cv_FCnet(y= y, x= x,
                 alpha= alpha, lambda= lambda,
                 nfolds= nfolds, rep_cv = rep_cv,
                 type.measure = type.measure,
                 intercept= intercept,
                 standardize= standardize,
                 ...)

  alpha= pars[["alpha"]]

  lambda= pars[["lambda"]]

  #fit with the final parameters
  fit= glmnet(x, y,
              alpha= alpha,
              lambda= lambda,
              standardize= standardize,
              intercept= intercept,
              ...)

  #wrap up results
  res= list(fit= fit,
            alpha= alpha,
            lambda= lambda,
            y= y)

  if(fit_stats){

    #avoid coeffs are dropped in sparse matrices
    #set lasso-eliminated coefficients to 0

    # #all vars
    # all_rn= rownames(coef(fit))
    #
    # #rownames with values
    # rnwv= coef(fit)[,1] %in% coef(fit)@x
    #
    # cf= rep(0, length= length(all_rn))
    #
    # cf[rnwv]= coef(fit)@x

    # #so by now we have all coefficients
    # #regardless of if they are dropped (e.g. LASSO)
    # coeffs= data.frame(Feature= seq(1, length(all_rn)-1, 1),
    #                    Coefficient= cf[-1])
    #

    #not sure why the hell did I need the chunk avove while I can simply do:

    cname= rownames(coef(fit))
    cf= coef(fit)[,1]

    coeffs= data.frame(Feature= cname,
                       Coefficient= cf)

    res[["coeffs"]]= coeffs

    #predict behavioral scores based on model
    p= predict(fit,
               s= lambda,
               newx= as.matrix(x))

    res[["predicted"]]= p

    gfstats= evalFCnet(true = y,
                       predicted = p)

    res[["fit_stats"]]= gfstats


  }

  return(res)

}
