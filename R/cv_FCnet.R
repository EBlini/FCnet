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
#' The crossvalidated lambda and alpha parameters are returned. For the model,
#' use the `FCnetLOO()` function to avoid overfitting.
#' Differently from `glmnet::cv.glmnet()`, here the
#' mean absolute error is minimized by default; you can change this parameter
#' directly in the function or by running `optionsFCnet(optionsFCnet= "mse")`.

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
#' times, and the consensus (by default: median) crossvalidated alpha and lambda values across all rep_cv
#' are returned. Useful in order to decrease randomness when kfold crossvalidation
#' is required. Setting this parameter too high though will result in a mere
#' approximation of the LOO.
#' @param cv.type.measure The measure to minimize in crossvalidation inner loops.
#' Differently from `glmnetUtils::cva.glmnet()` the deafult is the mean absolute error.
#' @param intercept whether to fit (TRUE) or not (FALSE) an intercept to the model.
#' @param standardize Whether x must be standardized internally to glmnet.
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

  if(length(alpha)== 1 & length(lambda)==1){
    stop("Scalar values for alpha and lambda have been supplied: a vector is needed.")
  }

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

    #return fit to make predictions quicker/avoid too much nesting?
    min_error= which.min(as.numeric(pars[rownames(pars)== "error"]))
    best= as.numeric(pars[rownames(pars)== "best"])[min_error]

    fit= cv_ridge[[min_error]]$modlist[[best]]

    #now hyperparameters
    lambda= as.numeric(pars[rownames(pars)== optionsFCnet("whichLambda")])

    #if more folds/repetitions, take the median as consensus
    if(rep_cv>1){

      lambda= optionsFCnet("consensus_function")(lambda)

    }

    alpha= as.numeric(pars[rownames(pars)== "alpha"])
    if(rep_cv>1){

      alpha= optionsFCnet("consensus_function")(alpha)

    }

    #return coefficients
    cname= rownames(coef(fit))
    cf= coef(fit)[,1]

    coeffs= data.frame(Feature= cname,
                       Coefficient= cf)

    rownames(coeffs)= NULL
    #useless at this point
    # #predict behavioral scores based on model
    # p= predict(fit,
    #            s= lambda,
    #            newx= as.matrix(x)
    #            )
    #
    # #fit stats - maybe unnecessary
    # gfstats= evalFCnet(true = y,
    #                    predicted = p)

    #return best parameters
    bp= list(alpha= alpha,
             lambda= lambda,
             fit= fit,
             y= y,
             coeffs= coeffs#,
             #predicted= p,
             #fit_stats= gfstats
             )


    return(bp)

  } #end if | vectors>1

}
