
#' Reduce the dimensionality of Functional Connectivity matrices or volumes
#'
#' This function accepts a list of lists containing squared FC matrices
#' or three-dimensional arrays (brain volumes) - such as one created by `loadFC()` - and apply feature reduction techniques.
#' Available techniques are listed under the `method` parameter and currently
#' include Principal Component Analysis (PCA) and Independent Component
#' Analysis (ICA). Both are applied on the upper triangular part of the matrices,
#' if matrices are provided, or to the whole array: if many brain volumes
#' are provided, consider masking them beforehand to reduce memory consumption.
#' PCA uses base R `prcomp()` while ICA uses `ica::icafast()`,
#' thus the `ica` package must be installed. The user can supply the desired
#' number of components to retain; else, the number of components needed to
#' explain at least a given proportion of the variance of FC matrices will be returned under
#' the `Weights` slot - for ICA, this is based on the PCA analysis. The default
#' is 95% of the varaince.


#' @param FCmatrices a list of lists including (squared) FC matrices
#' or three-dimensional arrays, such as
#' one provided by `loadFC()`.
#' @param method One of "PCA" or "ICA", defining the method for the desired
#' feature reduction technique.
#' @param Ncomp Number of components to retain. The default (NULL) automatically
#' retains 95% of the explained variance. If Ncomp== "all" returns all
#' the components. If Ncomp <1 this is interpreted as if the user wishes to retain
#' a given proportion of variance (e.g. 0.6). If Ncomp is a vector, an additional
#' slot is returned: a list of lists offering solutions with many components
#' (redundant for PCA but possibly useful for crossvalidating the number
#' of features in the case of ICA). In that last case the Weight slot always returns the model
#' accounting for at least 95% of variance.
#' @param parallel If TRUE uses `future.apply::future_lapply()`:
#' `future.apply` must be installed, your machine should have multiple cores
#' available for use, and threads should be defined explicitly
#' beforehand by the user (e.g. by calling `plan(multisession)`). Provided
#' these conditions are met, this can sensibly speed up the process - especially
#' for ICA estimations - when solutions with different components are asked
#' (i.e. Ncomp is a vector).
#' @param explore If TRUE only returns the proportion of explained variance
#' of all components, possibly useful to inform analyses.
#' @param ... More commands passed to `ica::icafast()`.

#' @return Weights is a matrix of reduced features. Loadings contains either the
#' eigenvectors or the estimated sources mapping the reduced features onto the
#' original space. SummaryPCA reports the amount of explained variance of the features
#' with respect to variability in the original matrices. MeanFC is the mean
#' FC matrix or volume. CrossWeights (optional) is a list of lists reporting weights and
#' loadings for a range of provided Ncomp. Dim is a dataframe with
#' information about the original dimensions.

#' @export


reduce_featuresFC= function(FCmatrices,
                            method= c("PCA", "ICA"),
                            Ncomp= NULL,
                            parallel= F,
                            explore= F,
                            ...){

  method= match.arg(method)

  #is Ncomp a vector?
  if(! is.null(Ncomp) & length(Ncomp)>1){
    Ncomp_vector= Ncomp
    Ncomp= NULL
  } else {

    Ncomp_vector= NULL

  }

  #First collapse the upper triangles -
  #only if matrices are provided
  #else work on the array
  if(class(FCmatrices[[1]])[1] %in% c("matrix", "data.frame")){

    original_dimensions= data.frame(type= "matrix",
                                    d1= nrow(FCmatrices[[1]]),
                                    d2= ncol(FCmatrices[[1]]))

  } else {

    original_dimensions= data.frame(type= "volume",
                                    d1= dim(FCmatrices[[1]])[1],
                                    d2= dim(FCmatrices[[1]])[2],
                                    d3= dim(FCmatrices[[1]])[3])
  }

  #if matrices
  if(original_dimensions[1,1]== "matrix"){

    if(parallel){

      FC= future.apply::future_lapply(FCmatrices, function(x){

        x[upper.tri(x, diag = F)]

      }, future.seed = T)

    } else {

      FC= lapply(FCmatrices, function(x){

        x[upper.tri(x, diag = F)]

      })
    }
  } else {#end if matrix}

    FC= lapply(FCmatrices, function(x){

      c(x)

    })

  }

  #then merge in columns
  FC= do.call(rbind, FC)


  #first PCA - which is needed anyway
  PCA= prcomp(FC)

  summaryPCA= summary(PCA)$importance

  Ncomp95= as.numeric(which(summaryPCA[3,]>0.95)[1])

  if(explore){

    print(paste("Components needed to explain 95% of the variance:",
                Ncomp95))

    return(summaryPCA)

  } else {

    #Ncomp changes accordingly
    if(is.null(Ncomp)){

      Ncomp= Ncomp95

    }

    if(Ncomp== "all"){

      Ncomp= length(summaryPCA[3,])

    }

    if (Ncomp<1){

      Ncomp= as.numeric(which(summaryPCA[3,]>Ncomp)[1])

    }

    if (Ncomp>length(summaryPCA[3,])){

      warning("You asked for too many components!")

    }

    if(method== "PCA"){

      #you already have the PCA object
      Loadings= PCA$rotation[,1:Ncomp]
      Weights= predict(PCA,
                       newdata= FC)[,1:Ncomp]

      if(!is.null(Ncomp_vector)){

        if(parallel){

          CrossWeights= future.apply::future_lapply(Ncomp_vector, function(x){

            L= PCA$rotation[,1:x]
            W= predict(PCA,
                       newdata= FC)[,1:x]
            return(list(Weights= W,
                        Loadings= L))

          }, future.seed = T)

        } else {

          CrossWeights= lapply(Ncomp_vector, function(x){

            L= PCA$rotation[,1:x]
            W= predict(PCA,
                       newdata= FC)[,1:x]

            return(list(Weights= W,
                        Loadings= L))

        })

        }

        names(CrossWeights)= paste0("PCA_", Ncomp_vector)

      } #end if vector

    } else {

      #method: "ICA"
      #preprocess
      ppFC= t(FC)
      ppFC= ppFC - rowMeans(ppFC)

      #if a vector was provided
      if(!is.null(Ncomp_vector)){

        #get CrossWeights
        if(parallel){

          CrossWeights= future.apply::future_lapply(Ncomp_vector, function(x){

            res= ica::icafast(ppFC, nc = x, ...)

            L= res$S
            W= res$M

            colnames(L)= paste0("ICA_", 1:x)
            colnames(W)= paste0("ICA_", 1:x)

            return(list(Weights= W,
                        Loadings= L))

          }, future.seed = T)

        } else {

          CrossWeights= lapply(Ncomp_vector, function(x){

            res= ica::icafast(ppFC, nc = x, ...)

            L= res$S
            W= res$M

            colnames(L)= paste0("ICA_", 1:x)
            colnames(W)= paste0("ICA_", 1:x)

            return(list(Weights= W,
                        Loadings= L))

          })


        } #end if parallel

        #now, tricky to match Ncomp: is it within the original vector?
        if(Ncomp %in% Ncomp_vector){
          Weights= CrossWeights[[which(Ncomp_vector== Ncomp)]]$Weights
          Loadings= CrossWeights[[which(Ncomp_vector== Ncomp)]]$Loadings
        } else {

          Weights= NULL
          Loadings=NULL

        }

        names(CrossWeights)= paste0("ICA_", Ncomp_vector)

      } # end if vector

      if(is.null(Ncomp_vector) | !Ncomp %in% Ncomp_vector){

        #do it
        ICA= ica::icafast(ppFC, nc = Ncomp, ...)

        Loadings= ICA$S
        Weights= ICA$M
        colnames(Loadings)= paste0("ICA_", 1:Ncomp)
        colnames(Weights)= paste0("ICA_", 1:Ncomp)



      }


    } #end method

    #merge and return
    res= list(Weights= Weights,
              Loadings= Loadings,
              summaryPCA= summaryPCA,
              dim= original_dimensions)

    #also return mean FC matrix or mean volume

    if(original_dimensions[1,1]== "volume"){

      meanFC= apply(FC, 2, mean)
      meanFC= array(meanFC, dim=c(original_dimensions[1,2],
                                  original_dimensions[1,3],
                                  original_dimensions[1,4]))

    } else { #if not a volume

    if(parallel){

      meanFC= Reduce("+",
                     future.apply::future_lapply(FCmatrices, function(x) {

                       as.matrix(x)/ length(FCmatrices)

                     }, future.seed = T))

    } else {

      meanFC= Reduce("+",
                   lapply(FCmatrices, function(x) {

                     as.matrix(x)/ length(FCmatrices)

                     }))
    }}

    res[["MeanFC"]]= meanFC

    #merge crossweights if vector was provided

    if(!is.null(Ncomp_vector)){

      res[["CrossWeights"]]= CrossWeights

    }

    #return list (of lists)
    return(res)

  } #end if explore

}
