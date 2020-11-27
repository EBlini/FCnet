#' CrossValidate FCnet model
#'
#' This function is a wrapper around `glmnetUtils::cva.glmnet()`, which is
#' itself built around `glmnet::cv.glmnet()`. For extended documentation,
#' the readers are encouraged to consult the respective pages.
#' `cv_FCnet()` requires two objects at minimum: `y` is a vector or data.frame
#' with exactly one column, corresponding to the (behavioral) score to predict; `x`
#' is a data.frame or a list of lists with an entry named "Weights",
#' which includes the independent variables. `x` can be - and is meant to be -
#' one object created by `reduce_featuresFC()`, but this is not strictly necessary.
#' The crossvalidated lambda and alpha parameters only are returned. For the model,
#' use the `FCnet()` function. Differently from `glmnet::cv.glmnet()`, here the
#' mean absolute error is minimized by default; you can change this parameter
#' directly in the function or by running `optionsFCnet(optionsFCnet= "mse")`.
#' Another difference in the defaults setting is that x is not standardized
#' in the call because we assume FC data to be already normalized.
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
#' @param standardize Whether x must be standardized. Differently from
#' `glmnet::glmnet()` the default is FALSE as we assume predictors are already either
#' summarised with PCA or ICA (and therefore scaled) or drawn from normalized FC matrices.
#' @param ... Other parameters passed to `glmnetUtils::cva.glmnet()`.
#'
#' @return The crossvalidated alpha and lambda parameters with the associated error.

#' @export



cv_FCnet= function(y, #dependent variable, typically behavior
                  x, #independent variables, typically neural measures
                  alpha= seq(0, 1, by= 0.1),
                  lambda= rev(10^seq(-5, 5, length.out = 200)),
                  nfolds= nrow(x),
                  rep_cv= 1,
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


  #to add here: what to do with missing values?


  #crossvalidate alpha and/or lambda, but only if vectors are supplied
  if (length(alpha)>1 | length(lambda)>1){

    #crossvalidate n times, take the median of the recovered parameters
    cv_ridge= lapply(1:rep_cv, function(time){

      cva= cva.glmnet(x= x, y= y,
                      alpha = alpha,
                      lambda = lambda,
                      nfolds= nfolds,
                      grouped= ifelse(nfolds== nrow(x), F, T),
                      type.measure = type.measure,
                      intercept= intercept,
                      standardize= standardize,
                      ...
                      )
      return(cva)
    })

    pars= sapply(cv_ridge, function(p)get_CVparsFCnet(p))

    lambda= as.numeric(pars[rownames(pars)== optionsFCnet("whichLambda")])

    alpha= as.numeric(pars[rownames(pars)== "alpha"])

    #return best parameters
    bp= list(alpha= alpha,
             lambda= lambda)

    #return fit to make predictions quicker/avoid nesting?
    best= as.numeric(pars[rownames(pars)== "best"])

    fit= cv_ridge[[1]]$modlist[[best]]

    bp[["fit"]]= fit


    return(bp)

  } #end if | vectors>1

}
