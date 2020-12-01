#' @title A set of functions for the simulation of Functional Connectivity matrices
#' @name simulateFCnet
#'
#' @description Provides a set of handy but crude functions in order to simulate
#' FC matrices starting from a model. Functions for simulating networks perturbations
#' are provided as well with some strict assumptions.
#'
#'
#' @param mat A square matrix used as template
#' @param matrices List of lists of FC matrices such as the one
#' created by `simulateMat()`.
#' @param Nmat Number of matrices to simulate from the template
#' @param mat_variability SD of the gaussian noise used to simulate
#' different matrices
#' @param y the behavioral score according to which networks will be
#' biased or perturbed. Recommended to be scaled.
#' @param network1 Indices of the first network to perturb
#' @param network2 Indices of the second network to perturb. If indices are
#' the same as those of network1 signal will be injected within one network, else
#' signal will be injected into the interaction between networks.
#' @param bias_multiplier The network will be perturbed according to this parameter.
#' At state, bias multiplier will be added or subtracted according to the
#' individual performance in y.
#' @param bias_variability Add a randomly distributed value with mean 0 and sd
#' bias_variability to the previous parameters such that participants may present
#' variability of predictivity within the perturbed network.
#' @param noise_multiplier Adds further gaussian noise within a network, with mean
#' 0 and sd= noise_multiplier, according to the individual performance in y. Worst
#' performances are injected with the entire noise_multiplier variability. Best
#' performances are injected with a minor variability.
#'
#'
#' @return A squared matrix or a list of matrices which have been polished after a simulation
#' with, for example, random variability or systematic drift.
#'

#'
#' @description **polishMat()** Polishes a matrix which has been subject to weird
#' shenanigans. Values are re-confined to the range c(-1, 1) and the
#' diagonal is restored.
#' @rdname simulateFCnet
#' @export

polishMat= function(mat){

  #copy upper triangle to lower triangle
  mat[lower.tri(mat)]=  t(mat)[lower.tri(mat)]

  #set to 1 correlations larger than 1
  mat[mat>1]= 1
  mat[mat<-1]= -1

  #diagonal back to 1
  diag(mat)= 1

  return(mat)

}

#'
#' @rdname simulateFCnet
#' @description **simulateMat()** Simulates Nmat matrices from a template matrix by
#' adding gaussian noise with sd = mat_variability.
#' @export
#this function creates several matrices by adding gaussian noise
#with sd "variability" to a reference matrix
simulateMat= function(mat,
                      Nmat,
                      mat_variability){

  sim= lapply(1:Nmat, function(s){

    ms= mat + rnorm(nrow(mat)*ncol(mat),
                    mean = 0, sd= mat_variability)

    ms= polishMat(ms)

  })

  return(sim)

}

#'
#' @rdname simulateFCnet
#' @description **biasMat()** Adds to a specific network of a list of matrices a systematic drift
#' which is a function of a given behavioral score.

#' @export
#this function biases a given network as a function
#of a behavioral score y
biasMat= function(matrices, y,
                  network1, network2,
                  bias_multiplier, bias_variability){

  sim= lapply(1:length(matrices), function(x){

    mat= matrices[[x]]

    #the score will be proportional to participants' performance
    #though with a random component included
    score= (y[x] / max(abs(y))) *
      (bias_multiplier + rnorm(1, 0, bias_variability))

    dim_mat= dim(mat[network1, network2])

    mat[network1, network2] = mat[network1, network2] +
      rep(score, dim_mat[1]*dim_mat[2])
    mat[network2, network1] = mat[network2, network1] +
      rep(score, dim_mat[1]*dim_mat[2])

    mat= polishMat(mat)

    return(mat)

  })

  return(sim)

}


#' @rdname simulateFCnet
#' @description **noiseMat()** Adds to a specific network of a list of matrices a randomly-distributed
#'  noise which is a function of a given behavioral score.

#' @export
#this function adds gaussian noise to a network
#amount of noise is proportional to y
noiseMat= function(matrices, y,
                   network1, network2,
                   noise_multiplier){

  lapply(1:length(matrices), function(x){

    mat= matrices[[x]]

    #variabilitry of the score (must be >0) function of y
    score= (y[x] + abs(min(y)))/max(abs(y)) + 0.001

    dim_mat= dim(mat[network1, network2])

    mat[network1, network2] = mat[network1, network2] +
      rnorm(dim_mat[1]*dim_mat[2], mean= 0, sd = score*noise_multiplier)

    mat[network2, network1] = mat[network2, network1] +
      rnorm(dim_mat[1]*dim_mat[2], mean= 0, sd = score*noise_multiplier)

    mat= polishMat(mat)

    return(mat)

  })

}

