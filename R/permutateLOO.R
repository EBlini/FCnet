#' Permutation test for FCnetLOO
#'
#' This function is a wrapper around `FCnet::FCnetLOO()` which executes nperm
#' times LOO robust regression pipelines on randomly permutated y scores.
#' A vector of R2 obtained from permutated (null) models is returned. If asked,
#' a data.frame (possibly huge) of the coefficients for the null models is also
#' returned. A summary data.frame describes the distribution of the
#' permutated models.
#'
#' @param y The dependent variable, typically behavioral scores to predict.
#' This can be a vector or a single data.frame column.
#' @param x The independent variables, typically neural measures that have
#' been already summarised through data reduction techniques
#' (e.g. ICA, PCA): an object created by `reduce_featuresFC()` will do. If such
#' an object is passed to this function, the "Weights" slot is taken as x.
#' Another kind of list can be passed to this function: in this case the function needs an
#' entry named "Weights". Otherwise, a data.frame can be passed to x.
#' @param alpha Value(s) that bias the elastic net toward ridge regression
#' (alpha== 0) or LASSO regression (alpha== 1). If a vector of alpha values
#' is supplied, the value is optimized through crossvalidation.
#' It defaults to a vector ranging from 0 to 1 with steps of 0.1.
#' @param lambda Regularization parameter for the regression,
#' see `glmnet::glmnet()`. Lambda must be a vector with length>1.
#' When a vector of lambda values is supplied, the value of lambda
#' is optimized through internal crossvalidation. It defaults to a vector
#' ranging from 10^-5 to 10^5 with 200 values in logarithmic steps.
#'
#' @param parallelLOO If TRUE - recommended, but not the default - uses
#' `future.apply::future_lapply()` for the outer loops: `future.apply` must be
#' installed, the machine should have multiple cores available for use,
#' and threads should be defined explicitly beforehand by the user
#' (e.g. by calling `plan(multisession)`).
#' @param scale_y Whether y should be scaled prior to fit. Default, TRUE, scales
#' and center y with `scale()`.
#'
#' @param nperm The number of permutations for the null models. Default is 100.
#'
#' @param model_R2 Optional. If this entry is left NULL, the original model is
#' fitted again. Either an object created by `FCnet::FCnetLOO()` or a
#' precise value of R2 can be supplied, as to avoid unnecessary computations.
#'
#' @param return_coeffs Optional: whether coefficients for the null models should
#' be returned as well. This may interesting should inferential statistics be
#' envisaged for single coefficients. The returned datya.frame, on the other
#' hand, may be quite large.

#' @param cv.type.measure The measure to minimize in crossvalidation inner loops.
#' Differently from `glmnetUtils::cva.glmnet()` the default is the mean absolute error.
#' @param intercept whether to fit (TRUE) or not (FALSE) an intercept to the model.
#' @param standardize Whether x must be standardized internally to glmnet.
#' @param ... Other parameters passed to `glmnetUtils::cva.glmnet()` or
#' `glmnet::glmnet()`.
#'
#' @return A list including the R2 for the permutated models and a summary data.frame. Optionally,
#' a data.frame including the permutated coefficients.
#'
#' @export

permutateLOO= function(y,
                       x,
                       alpha= seq(0, 1, by= 0.1),
                       lambda= rev(10^seq(-5, 5, length.out = 200)),
                       parallelLOO= F,
                       scale_y= T,
                       nperm= 100,
                       model_R2= NULL,
                       return_coeffs= F,
                       type.measure= optionsFCnet("cv.type.measure"),
                       intercept= optionsFCnet("intercept"),
                       standardize= optionsFCnet("standardize"),
                       ...){

  #ensure you are working with matrices
  y= data.matrix(y)

  #scaling if requested, then set to False in inner call
  if(scale_y){y= scale(y)}

  #original model here
  if(is.null(model_R2)){

    original_R2= FCnetLOO(y= y,
                          x= x,
                          alpha = alpha,
                          lambda= lambda,
                          parallelLOO= parallelLOO,
                          scale_y= F,
                          type.measure= type.measure,
                          intercept= intercept,
                          standardize= standardize,
                          ...)$R2

  } else {

    if(class(model_R2)[1]== "list")(original_R2= model_R2$R2)

  }

  if(is.na(original_R2))(original_R2= 0)
  original_R2= as.numeric(original_R2)

  if(parallelLOO){

    perms= future.apply::future_lapply(1:nperm, function(p){

      y= y[sample(1:length(y), size = length(y), replace = F)]

      res= FCnetLOO(y= y,
                    x= x,
                    alpha = alpha,
                    lambda= lambda,
                    parallelLOO= F,
                    scale_y= F,
                    type.measure= type.measure,
                    intercept= intercept,
                    standardize= standardize,
                    ...
                    )

      R2= data.frame(nPerm= p, R2= ifelse(is.na(res$R2), 0, res$R2))

      res_inner= list(R2=R2)

      if(return_coeffs){

        Coeffs= res$coeffs
        Coeffs$nPerm= p
        res_inner[["Coeffs"]]= Coeffs

      }

      return(res_inner)

    }, future.seed= T)} else {

      perms= lapply(1:nperm, function(p){

        y= y[sample(1:length(y), size = length(y), replace = F)]

        res= FCnetLOO(y= y,
                      x= x,
                      alpha = alpha,
                      lambda= lambda,
                      parallelLOO= F,
                      scale_y= F,
                      type.measure= type.measure,
                      intercept= intercept,
                      standardize= standardize,
                      ...
                      )

        R2= data.frame(nPerm= p, R2= res$R2)

        res_inner= list(R2=R2)

        if(return_coeffs){

          Coeffs= res$coeffs
          Coeffs$nPerm= p
          res_inner[["Coeffs"]]= Coeffs

        }

        return(res_inner)

      })

    } #end if parallelLOO


  R2= lapply(perms, function(x)x$R2)
  R2= do.call(rbind, R2)

  res= list(R2= R2)

  if(return_coeffs){

    Coeffs= lapply(perms, function(x)x$Coeffs)
    Coeffs= do.call(rbind, Coeffs)
    res[["Coeffs"]]= Coeffs
  }

  #summary
  summary_perms= data.frame(Model_R2= original_R2,
                            Permutations= nperm,
                            P_value= sum(R2$R2>original_R2)/length(R2$R2),
                            MeanR2= mean(R2$R2),
                            SD_R2= sd(R2$R2),
                            Upper_95= mean(R2$R2) + 1.96*sd(R2$R2),
                            Quantile_50= quantile(R2$R2, 0.5),
                            Quantile_70= quantile(R2$R2, 0.7),
                            Quantile_90= quantile(R2$R2, 0.9),
                            Quantile_95= quantile(R2$R2, 0.95))

  rownames(summary_perms)= NULL

  res[["Summary"]]= summary_perms

  return(res)

}
