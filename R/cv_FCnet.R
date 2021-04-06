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
#' @param parallelLOO If TRUE - recommended, but not the default - uses
#' `future.apply::future_lapply()` for the outer loops: `future.apply` must be
#' installed, the machine should have multiple cores available for use,
#' and threads should be defined explicitly beforehand by the user
#' (e.g. by calling `plan(multisession)`).
#' @param cv_Ncomp Whether to crossvalidate the number of components or not.
#' It defaults to NULL, but a vector can be supplied specifing the number (range) of
#' components to test in the inner loops.
#' @param cv_Ncomp_method Whether the number of components to optimize means
#' components are ordered (e.g. according to the explained variance of neuroimaging
#' data) or - somehow experimental - whether to use the N best components
#' ranked according to their relationship (pearson's R) with y.
#' @param cv.type.measure The measure to minimize in crossvalidation inner loops.
#' Differently from `glmnetUtils::cva.glmnet()` the deafult is the mean absolute error.
#' @param intercept whether to fit (TRUE) or not (FALSE) an intercept to the model.
#' @param standardize Whether x must be standardized internally to glmnet.
#' @param thresh Threshold for glmnet to stop converging to the solution.
#' @param ... Other parameters passed to `glmnetUtils::cva.glmnet()`.
#'
#' @return The crossvalidated alpha and lambda parameters with the associated error.

#' @export

cv_FCnet= function(y, #dependent variable, typically behavior
                  x, #independent variables, typically neural measures
                  alpha= seq(0, 1, by= 0.1),
                  lambda= rev(10^seq(-5, 5, length.out = 200)),
                  parallelLOO= F,
                  cv_Ncomp= NULL,
                  cv_Ncomp_method= c("order", "R"),
                  type.measure= optionsFCnet("cv.type.measure"),
                  intercept= optionsFCnet("intercept"),
                  standardize= optionsFCnet("standardize"),
                  thresh= optionsFCnet("thresh"),
                  ...){


  cv_Ncomp_method= match.arg(cv_Ncomp_method)


  if(length(lambda)==1){
    stop("Scalar values for lambda have been supplied: a vector is needed.")
  }

  #ensure you are working with matrices
  y= data.matrix(y)

  if(class(x)[1]== "list"){

    x= data.matrix(x$Weights)

  } else {

    if(class(x)[1]== "data.frame"){x= data.matrix(x)}

  }

  #nfold is set to loo
  nfolds= nrow(x)

  #initialize a vector of zeros to store coefficients
  #this is in case cv of the number of components is required
  cname= c(("Intercept"), colnames(x))
  zeros= rep(0, ncol(x) + 1) #+1 is the intercept

  #if cv_Ncomp is null, all components are tested
  if(is.null(cv_Ncomp))(cv_Ncomp= ncol(x))

  if(cv_Ncomp_method== "order" | length(cv_Ncomp)== 1){

    test_c= lapply(cv_Ncomp, function(x)1:x)

  }

  if(cv_Ncomp_method== "R") {

    r= apply(x, 2, function(z){cor(z, y)})

    rank= rank(-abs(r), ties.method = "max")
    all_x= 1:ncol(x)

    test_c= lapply(cv_Ncomp, function(n){

      all_x[rank<= n]

    })
  }

  #to add here: what to do with missing values?


  #crossvalidate alpha and/or lambda, but only if vectors are supplied
  if (length(alpha)>1 | length(lambda)>1 | length(cv_Ncomp)>1){

    if(parallelLOO== T){

      Ncomp_ridge= future.apply::future_lapply(test_c, function(NC){

        cva= cva.glmnet(x= x[, NC], y= y,
                        alpha = alpha,
                        lambda = lambda,
                        nfolds= nfolds,
                        grouped= ifelse(nfolds== nrow(x), F, T),
                        type.measure = type.measure,
                        intercept= intercept,
                        standardize= standardize,
                        thresh= thresh,
                        ...
        )
        return(cva)

      }, future.seed= T)

    } else {

      Ncomp_ridge= lapply(test_c, function(NC){

        cva= cva.glmnet(x= x[, NC], y= y,
                        alpha = alpha,
                        lambda = lambda,
                        nfolds= nfolds,
                        grouped= ifelse(nfolds== nrow(x), F, T),
                        type.measure = type.measure,
                        intercept= intercept,
                        standardize= standardize,
                        thresh= thresh,
                        ...
        )
        return(cva)

      })
    }




    pars= sapply(Ncomp_ridge, function(p)get_CVparsFCnet(p))

    #the minimum error observed among the external loop, i.e. ncomp
    min_error= which.min(as.numeric(pars[rownames(pars)== "error"]))

    #the best model (alpha, etc) where the minimum error is observed
    best= as.numeric(pars[rownames(pars)== "best"])[min_error]

    fit= Ncomp_ridge[[min_error]]$modlist[[best]]


    #now hyperparameters of the best model,
    #i.e. where the minimum error was observed
    lambda= as.numeric(pars[rownames(pars)== optionsFCnet("whichLambda")])[min_error]
    alpha= as.numeric(pars[rownames(pars)== "alpha"])[min_error]
    N_comp= as.numeric(pars[rownames(pars)== "N_comp"])[min_error]

    which_comp= test_c[[min_error]]

    #return coefficients - only for components that were actually tested
    #the rest is set to zero
    cf= coef(fit, s= lambda)[,1]
    zeros[c(1, which_comp+1)] = cf
    cf= zeros

    coeffs= data.frame(Feature= cname,
                       Coefficient= cf)

    rownames(coeffs)= NULL


    #return best parameters
    bp= list(alpha= alpha,
             lambda= lambda,
             N_comp= which_comp,
             fit= fit,
             y= y,
             coeffs= coeffs
             )


    return(bp)

  } #end if | vectors>1

}
