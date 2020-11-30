
#' Back-projection of regression coefficients onto the original space
#'
#' This function projects - based on Loadings as provided,
#' for example, by `reduce_features_FC()` - regression coefficients onto
#' the original space (i.e. a square matrix). Coefficients can either be
#' passed as a list (as created, for example, by `FCnetLOO()`), as a data.frame,
#' or as a user-defined vector. If a list or a data.frame is provided, this function
#' finds consensus trough the function specified by `optionsFCnet("consensus_function")`,
#' which by default is the median. Moreover, if a single value for the coefficients
#' is provided, then the function understands that the user wants to retrieve
#' the edges' loadings for a specific coefficients, i.e. (possibly) a network.
#' The function returns a square matrix whose entries are the linear combination
#' of the coefficients and the edge's loadings. Optional thresholding can be passed
#' to the threshold parameter: an integer value must be specified in order to
#' retrieve the k largest coefficients, in absolute value, only.
#'
#' @param coeffs Coefficients can either be
#' passed as a list (as created, for example, by `FCnetLOO()`), as a data.frame,
#' or as a user-defined vector. If a list or a data.frame is provided, this function
#' finds consensus trough the function specified by `optionsFCnet("consensus_function")`,
#' which by default is the median. Moreover, if a single value for the coefficients
#' is provided, then the function understands that the user wants to retrieve
#' the edges' loadings for a specific coefficients, i.e. (possibly) a network.
#' @param reduce_features_object An object created by `reduce_features_FC()`. If present,
#' Loadings may not be provided, and is overwritten.
#' @param Loadings Not necessary if reduce_features_object is passed to the function.
#' Else, this is a dataframe containing all edges' relationships with the coefficients.
#'@param threshold Optional. Prune the back-projection matrix by retaining only the
#'threshold larger entries (in absolute value).
#'
#' @return A square back-projection matrix with the original dimensions.
#'
#' @export
#'

backprojectFCnet= function(coeffs,
                           reduce_features_object,
                           Loadings= NULL,
                           threshold= NULL){


  #work with coefficients
  if(class(coeffs)== "list")(coeffs= coeffs$coeffs)
  if(class(coeffs)== "data.frame"){

    coeffs= tapply(coeffs$Coefficient,
                   coeffs$Feature,
                   optionsFCnet("consensus_function"))

    if("(Intercept)" %in% names(coeffs)){

      intercept= as.numeric(coeffs[names(coeffs) %in% "(Intercept)"])
      coeffs= as.numeric(coeffs[!names(coeffs) %in% "(Intercept)"])

    }
  }

  #Weights and Loadings
  if(!is.null(reduce_features_object)){
    Loadings= reduce_features_object$Loadings
  }

  #if coeff is a single number - depends on length wights components
  if(length(coeffs)== 1){
    c= rep(0, ncol(Loadings))
    c[coeffs]= 1
    coeffs= c
  }

  #warn if all coeffs are 0
  make_sense= sum(coeffs==0) != length(coeffs)

  if(!make_sense){warning("All coefficients are 0: model is not predictive.")}


  #do the job here

  #multiply coeffs per PCA loadings
  final_coeffs= t(t(Loadings) * coeffs)


  final_coeffs= rowSums(final_coeffs)

  # final_coeffs= Weights %*% t(final_coeffs)
  #
  # #find linear contribution for each row (connection)
  # final_coeffs= colSums(final_coeffs)


  my_backTransform= function(v){

    bt1= matrix(c(0),
                ceiling(sqrt(1.0+8.0*length(v))/2.0),
                ceiling(sqrt(1.0+8.0*length(v))/2.0))

    bt1[upper.tri(bt1, diag=FALSE)]= v

    bt1[lower.tri(bt1)]=  t(bt1)[lower.tri(bt1)]

    return(bt1)

  }

  #vector to matrix
  res= my_backTransform(final_coeffs)


  if(!is.null(threshold)){

    Nth_larger= sort(abs(final_coeffs),
                     decreasing = T)[threshold]

    final_coeffs_pruned= final_coeffs

    final_coeffs_pruned[abs(final_coeffs_pruned)<Nth_larger]= 0

    res= my_backTransform(final_coeffs_pruned)

  }

  return(res)


}
