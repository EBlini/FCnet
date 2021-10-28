
#from vector to matrix
vtom= function(v){

  bt1= matrix(c(0),
              ceiling(sqrt(1.0+8.0*length(v))/2.0),
              ceiling(sqrt(1.0+8.0*length(v))/2.0))

  bt1[upper.tri(bt1, diag=FALSE)]= v

  bt1[lower.tri(bt1)]=  t(bt1)[lower.tri(bt1)]

  return(bt1)

}

vtoarray= function(v, dim){

  bt1= array(v, dim)

  return(bt1)

}

#' Back-projection of regression coefficients onto the original space
#'
#' This function projects - based on Loadings as provided
#' by `reduce_features_FC()` - regression coefficients onto
#' the original space (i.e. a square matrix or a volume).
#' Coefficients can either be
#' passed as a list (as created, for example, by `FCnetLOO()`), as a data.frame,
#' or as a user-defined vector. If a single value for the coefficients
#' is provided, then the function understands that the user wants to retrieve
#' the edges' loadings for a specific coefficient, i.e. (possibly) a network.
#' The function returns a square matrix or a 3-D array whose entries are the linear combination
#' of the coefficients and the edge's loadings. Optional thresholding can be passed
#' to the threshold parameter: an integer value must be specified in order to
#' retrieve the k largest coefficients, in absolute value, only.
#'
#' @param coeffs Coefficients can either be
#' passed as a list (as created, for example, by `FCnetLOO()`), as a data.frame,
#' or as a user-defined vector. If you fetch a list or a data.frame,
#' then all the coefficients will be backprojected: you may want to select those that are
#' significant following a permutation test instead (see `select_coefficients()`).
#' If a single value for the coefficients
#' is provided, then the function understands that the user wants to retrieve
#' the edges' loadings for a specific coefficients, i.e. (possibly) a network.
#' @param reduce_features_object An object created by `reduce_features_FC()`. If present,
#' Loadings may not be provided, and is overwritten.
#' @param normthresh If TRUE first back-projected coefficients are normalized
#' within the -1 - 1 interval, and then only the values larger than 0.1 (in
#' absolute value) are retained. Default is TRUE.
#' @param threshold Optional. Prune the back-projection matrix by retaining only the
#'threshold larger entries (in absolute value) if normthresh is FALSE. This is disregarded
#'if normthresh is TRUE.
#' @return A square back-projection matrix or 3-D volume with
#' the original dimensions.
#'
#' @export
#'

backprojectFCnet= function(coeffs,
                           reduce_features_object,
                           normthresh= TRUE,
                           threshold= NULL){


  #work with coefficients
  if(class(coeffs)== "list")(coeffs= coeffs$coeffs)
  if(class(coeffs)== "data.frame"){

    #get coeffs names - useful for reordering them
    names= coeffs$Feature

    coeffs= coeffs$Coefficient

    # #reorder
    # coeffs= coeffs[names]


    if("Intercept" %in% names | "(Intercept)" %in% names){

      intercept= as.numeric(coeffs[names %in% c("(Intercept)", "Intercept")])
      coeffs= as.numeric(coeffs[!names %in% c("(Intercept)", "Intercept")])

    }
  }

  #Weights and Loadings
  Loadings= reduce_features_object$Loadings


  #if coeff is a single number -
  # or a subset of loadings
  #then depends on length weights components
  if(length(coeffs)< ncol(Loadings)){
    c= rep(0, ncol(Loadings))
    c[coeffs]= 1
    coeffs= c
  }

  #warn if all coeffs are 0
  make_sense= sum(coeffs==0) != length(coeffs)

  if(!make_sense){warning("All coefficients are 0: model is not predictive.")}


  #do the job here
  #multiply coeffs for PCA loadings
  final_coeffs= t(t(Loadings) * coeffs)


  #sum or mean?
  final_coeffs= rowSums(final_coeffs)

  #warn if both threshold pars are set
  if(!is.null(threshold) & normthresh== TRUE){

    warning("Parameter threshold disregarded because normthresh is TRUE")

    threshold= NULL

  }

  #vector to matrix
  if(is.null(threshold) & normthresh== FALSE){

    #no pruning
    if(reduce_features_object$dim[1,1]== "matrix"){

      res= vtom(final_coeffs)

    } else {

      res= vtoarray(final_coeffs, c(reduce_features_object$dim[1,2],
                                    reduce_features_object$dim[1,3],
                                    reduce_features_object$dim[1,4]))

    }

  } else if (!is.null(threshold) & normthresh== FALSE){

    #if pruning by threshold
    Nth_larger= sort(abs(final_coeffs),
                     decreasing = T)[threshold]

    final_coeffs_pruned= final_coeffs

    final_coeffs_pruned[abs(final_coeffs_pruned)<Nth_larger]= 0

    if(reduce_features_object$dim[1,1]== "matrix"){

      res= vtom(final_coeffs_pruned)

    } else {

      res= vtoarray(final_coeffs_pruned,
                    c(reduce_features_object$dim[1,2],
                      reduce_features_object$dim[1,3],
                      reduce_features_object$dim[1,4]))

    }


  } else {

    #if normtresh is TRUE
    larger= max(abs(final_coeffs))

    final_coeffs_pruned= final_coeffs/larger

    final_coeffs_pruned[abs(final_coeffs_pruned)< 0.1]= 0

    if(reduce_features_object$dim[1,1]== "matrix"){

      res= vtom(final_coeffs_pruned)

    } else {

      res= vtoarray(final_coeffs_pruned,
                    c(reduce_features_object$dim[1,2],
                      reduce_features_object$dim[1,3],
                      reduce_features_object$dim[1,4]))

    }

  }

  return(res)


}
