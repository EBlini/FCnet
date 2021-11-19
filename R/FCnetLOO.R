#' Leave-One-Out fit of elastic nets for Functional Connectivity data with nested crossvalidation
#'
#' This function is a wrapper around `glmnet::glmnet()` as called from
#' `FCnet::cv_FCnet()`. For extended documentation,
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
#' details of the crossvalidation procedures can be passed as `...` arguments to
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
#' is supplied, the value is optimized through nested crossvalidation.
#' It defaults to a vector ranging from 0 to 1 with steps of 0.1.
#' The crossvalidated alpha is returned.
#' @param lambda Regularization parameter for the regression,
#' see `glmnet::glmnet()`. Lambda must be a vector with length>1.
#' When a vector of lambda values is supplied, the value of lambda
#' is optimized through internal nested crossvalidation. It defaults to a vector
#' ranging from 10^-5 to 10^5 with 200 values in logarithmic steps.
#' The crossvalidated optimal lambda is returned.
#' @param cv_Ncomp Whether to crossvalidate the number of components or not.
#' It defaults to NULL, but a vector can be supplied specifing the number (range) of
#' components to test in the inner loops.
#' @param cv_Ncomp_method Whether the number of components to optimize means
#' components are ordered (e.g. according to the explained variance of neuroimaging
#' data) or - somehow experimental - whether to use the N best components
#' ranked according to their relationship (the coefficient of an
#' univariate (g)lm) with y.
#' @param parallelLOO If TRUE - recommended, but not the default - uses
#' `future.apply::future_lapply()` for the outer loops: `future.apply` must be
#' installed, the machine should have multiple cores available for use,
#' and threads should be defined explicitly beforehand by the user
#' (e.g. by calling `plan(multisession)`).
#' @param scale_y Whether y should be scaled prior to fit. Default, TRUE, scales and center y with
#' `scale()`.
#' @param scale_x Whether x should be scaled prior to fit. Default, TRUE, subtracts
#' the mean matrix value and divides each entry for the matrix variance.
#' Beware that this adds to `optionsFCnet("standardize")`.
#' @param family Defaults to "gaussian." Experimental support for "binomial" on the way.
#' @param cv.type.measure The measure to minimize in crossvalidation inner loops.
#' Differently from `glmnetUtils::cva.glmnet()` the default is the mean absolute error.
#' @param intercept whether to fit (TRUE) or not (FALSE) an intercept to the model.
#' @param standardize Whether x must be standardized.
#' @param thresh Threshold for glmnet to stop converging to the solution.
#' @param ... Other parameters passed to `glmnetUtils::cva.glmnet()` or
#' `glmnet::glmnet()`.
#'
#' @return Goodness of fit statistics for the outer loops, as well as LOO predictions.
#' Crossvalidated best alpha and lambda values for the inner loops as well as
#' all inner models' coefficients combined..
#'
#' @export


FCnetLOO= function(y,
                   x,
                   alpha= seq(0, 1, by= 0.1),
                   lambda= rev(10^seq(-5, 5, length.out = 200)),
                   cv_Ncomp= NULL,
                   cv_Ncomp_method= c("order", "R"),
                   parallelLOO= F,
                   scale_y= T,
                   scale_x= T,
                   family= optionsFCnet("family"),
                   type.measure= optionsFCnet("cv.type.measure"),
                   intercept= optionsFCnet("intercept"),
                   standardize= optionsFCnet("standardize"),
                   thresh= optionsFCnet("thresh"),
                   ...){

  cv_Ncomp_method= match.arg(cv_Ncomp_method)


  #scaling if requested
  if(scale_y){

    if(family== "binomial")(warning("You requested scaling of y but the family is 'binomial'..."))

    y= scale(y)

    }

  #ensure you are working with matrices
  y= data.matrix(y)


  #further reformatting of x
  if(class(x)[1]== "list"){

    x= data.matrix(x$Weights)

  } else {

    if(class(x)[1]== "data.frame"){x= data.matrix(x)}

  }

  #scaling if requested
  if(scale_x){x= (x - mean(x))/var(as.vector(x))}


  #main procedures here

  if (optionsFCnet("nested")== FALSE){

    #cv
    fit= cv_FCnet(y = y,
                  x = x,
                  alpha = alpha,
                  lambda = lambda,
                  parallelLOO= parallelLOO,
                  cv_Ncomp = cv_Ncomp,
                  cv_Ncomp_method = cv_Ncomp_method,
                  family= family,
                  type.measure= type.measure,
                  intercept= intercept,
                  standardize= standardize,
                  thresh= thresh,
                  ...
                  )

    p= fit$predictions

    rss <- sum((p - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss


    #metrics
    pars= evalFCnet(fit, family)


    #wrap-up info
    res= list(R2= rsq,
              Goodness_Fit= pars,
              predicted= p,
              alpha= fit$alpha,
              lambda= fit$lambda,
              N_comp= fit$N_comp,
              coeffs= fit$coeffs,
              y= y)

    if(family== "binomial"){

      res[["R2"]]= NULL

      res[["Accuracy"]]= pars$Accuracy

    }


   } else { # if  nested


  #indices for the LOO function below, passed within lapply or future_lapply
  lapply_over= 1:nrow(x)

  #main function here: parallel or not, nested
  loo_f= function(r){

    new_x= x[-r,]

    new_y= y[-r,]

    fit= cv_FCnet(y = new_y,
                  x = new_x,
                  alpha = alpha,
                  lambda = lambda,
                  cv_Ncomp = cv_Ncomp,
                  cv_Ncomp_method = cv_Ncomp_method,
                  parallelLOO = F,
                  family= family,
                  type.measure= type.measure,
                  intercept= intercept,
                  standardize= standardize,
                  thresh= thresh,
                  ...
    )

    p= predict(fit$fit,
               s= fit$lambda,
               newx= t(data.matrix((x[r, fit$N_comp]))),
               exact= TRUE)

    p= as.numeric(p)

    return(list(prediction= p,
                alpha= fit$alpha,
                lambda= fit$lambda,
                N_comp= length(fit$N_comp),
                coeffs= fit$coeffs))

  }



  if(parallelLOO== TRUE){

    loo= future.apply::future_lapply(lapply_over, function(r){

      loo_f(r)

      }, future.seed = T)

  } else {

    loo= lapply(lapply_over, function(r){

      loo_f(r)

      })

  } #end if parallelLOO


  ####################################################
  #### INNER #########################################
  ####################################################

  #extract and reshape all relevant information
  prediction_inner= sapply(loo, function(x)x[["prediction"]])

  alpha_inner= sapply(loo, function(x)x[["alpha"]])

  lambda_inner= sapply(loo, function(x)x[["lambda"]])

  N_comp_inner= sapply(loo, function(x)x[["N_comp"]])

  coeffs_inner= lapply(loo, function(x){data.frame(x$coeffs)})
  coeffs_inner= do.call(rbind, coeffs_inner)
  coeffs_inner$ID= rep(lapply_over,
                 each= nrow(coeffs_inner)/length(lapply_over))

  inner= list(predicted= prediction_inner,
              alpha= alpha_inner,
              lambda= lambda_inner,
              N_comp= N_comp_inner,
              coeffs= coeffs_inner)

  ####################################################
  #### OUTER #########################################
  ####################################################
  #find consensus
  consensus_alpha= optionsFCnet("consensus_function")(alpha_inner)
  consensus_lambda= optionsFCnet("consensus_function")(lambda_inner)
  consensus_N_comp= optionsFCnet("consensus_function")(N_comp_inner)

  #prepare
  fit= cv_FCnet(y = y,
                x = x,
                alpha = consensus_alpha,
                lambda = consensus_lambda + 0:1, #ensure it's a vector so cv.glmnet is happy
                cv_Ncomp = consensus_N_comp,
                cv_Ncomp_method = cv_Ncomp_method,
                family= family,
                type.measure= type.measure,
                intercept= intercept,
                standardize= standardize,
                thresh= thresh,
                ...
  )

  #recover number components + prepare coefficients
  cname= c(("Intercept"), colnames(x))
  zeros= rep(0, ncol(x) + 1) #+1 is the intercept

  #if cv_Ncomp is null, all components are tested
  if(is.null(cv_Ncomp))(cv_Ncomp= ncol(x))

  if(cv_Ncomp_method== "order" | length(cv_Ncomp)== 1){

    test_c= lapply(cv_Ncomp, function(x)1:x)

  }

  if(cv_Ncomp_method== "R") {

    if(family== "gaussian"){

      r= unlist(apply(x, 2, function(z){coef(lm(y~z))[2]}))

    } else {

      r= unlist(apply(x, 2, function(z){coef(glm(y~z, family= "binomial"))[2]}))

    }

    rank= rank(-abs(r), ties.method = "max")
    all_x= 1:ncol(x)

    test_c= lapply(cv_Ncomp, function(n){

      all_x[rank<= n]

    })
  }

  final_components= sapply(test_c, function(l)length(l)== consensus_N_comp)
  final_components= unlist(test_c[final_components])

  #CV predictions
  p= fit$predictions

  rss <- sum((p - y) ^ 2)  ## residual sum of squares
  tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss

  pars= evalFCnet(fit, family)

  #wrap-up info
  res= list(R2= rsq,
            Goodness_Fit= pars,
            predicted= p,
            alpha= consensus_alpha,
            lambda= consensus_lambda,
            N_comp= consensus_N_comp,
            coeffs= fit$coeffs,
            y= y,
            inner= inner)

  if(family== "binomial"){

    res[["R2"]]= NULL

    res[["Accuracy"]]= pars$Accuracy

  }


   } #end if nested

  return(res)

}
