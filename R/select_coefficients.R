#' Selection of coefficients based on the outcome of permutateLOO
#'
#' This function assesses the original coefficients against the ones obtained
#' from permutated samples as returned by `permutateLOO()` (with the parameter
#' `return_coeffs` set to `TRUE`). It returns a vector of weights as well as a
#' dichotomic mask reflecting which values pass a given threshold.
#'
#'
#' @param coeffs Coefficients can either be
#' passed as a list (as one created, for example, by `FCnetLOO()`), or as a data.frame (which must, however,
#' have the same structure as the one returned by `FCnetLOO()`).
#' This function finds consensus coefficients first trough the function specified
#' by `optionsFCnet("consensus_function")`,which by default is the median.
#' @param permutated_coeffs A list, such as an object created by
#' `permutateLOO()` (with the parameter `return_coeffs` set to `TRUE`),
#' or a data.frame (which must, however, have the same structure
#' as the one returned by `permutateLOO()`)
#' @param threshold Significance threshold against which evaluate the coefficients.
#' It is assumed to be two-tailed, and defaults to 0.05.
#'
#' @return A list with two entries: coeffs, i.e. the values of coefficients passing the
#' statistical threshold, and mask, reporting which coefficients are retained.
#'
#' @export

select_coefficients= function(coeffs,
                              permutated_coeffs,
                              threshold= 0.05){

  #work with coefficients
  if(class(coeffs)== "list")(coeffs= coeffs$coeffs)
  if(class(permutated_coeffs)== "list")(permutated_coeffs= permutated_coeffs$Coeffs)
  pc= permutated_coeffs #renamed for simplicity

  #get coeffs names - useful for reordering them
  names= coeffs$Feature

  #find consensus
  consensus= coeffs$Coefficient

  #find consensus for each permutation - not necessary but left for legacy
  perm_range= tapply(pc$Coefficient,
                     list(pc$Feature, pc$nPerm),
                     optionsFCnet("consensus_function"))

  #reorder
  perm_range= perm_range[names,]

  #find quantile with respect to permutated range
  mask= sapply(1:nrow(perm_range), function(x){

    min= min(perm_range[x,], consensus[x])

    q= ecdf(perm_range[x,]+abs(min))(consensus[x]+abs(min))

    m= ifelse(q < threshold/2 | q > (1 - threshold/2), 1, 0)

    #if they are all the same (zeros...)
    if(consensus[x]== 0 | sd(perm_range[x,])== 0)(m= 0)

    return(m)

  })

  coeffs= consensus
  coeffs[mask== 0]= 0

  #omit intercept
  if("Intercept" %in% names){
    coeffs= coeffs[-1]
    mask= mask[-1]
  }

  res= list(coeffs= coeffs,
            mask= mask)

  return(res)


}



