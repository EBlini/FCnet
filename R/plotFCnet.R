#' plot a FCnet model
#'
#' This function accepts a FCnet model - such as one retuned by
#' `FCnetLOO()` - and returns a plot for either the coefficients or the
#' model's prediction.
#'
#' @param model The FCnet model to depict.
#' @param output What to depict: model's predictions or models' coefficients. Coefficients
#' are assumed to arise from nested models, thus they will depicted with - more or less pronounced -
#' variability with the geometries chosen with the style parameter.
#' @param plot_labels if `TRUE` annotates the row names of x as predictions.
#' @param subset_coeffs To avoid padding, sometimes the user may want to subset the coefficients to report.
#' This can be done by passing a vector to this parameter. Note that the intercept counts as
#' a coefficient (of 0) even if omitted in the call.
#' @return A `ggplot2` object which can be customized further.

#' @export
#'

plotFCnet= function(model,
                    output= c("predictions", "coefficients"),
                    plot_labels= T,
                    subset_coeffs= NULL){


  output= match.arg(output)

  #commonTheme
  commonTheme= ggplot2::theme(
    text= ggplot2::element_text(size= 14, face= "bold"),
    axis.text= ggplot2::element_text(size= 14, face= "bold", color= "black"))


  if(output== "predictions"){

    ggDF= data.frame(ID= 1:length(model$predicted),
                     Score= as.numeric(model$y),
                     Predicted= model$predicted)

    #check if binomial
    if(length(table(ggDF$Score))==2)(binomial= T)else(binomial=F)

    lim= c(min(ggDF$Score, ggDF$Predicted),
           max(ggDF$Score, ggDF$Predicted))

    if(binomial) {

      ggDF$Score= as.factor(ggDF$Score)

      lim= c(-abs(max(ggDF$Predicted)),
                       abs(max(ggDF$Predicted)))

      }


    #plot here
    if(plot_labels){

      p= ggplot2::ggplot(ggDF,
                       ggplot2::aes(x= Score,
                                    y= Predicted,
                                    label= ID)) +
      ggplot2::geom_label()

    } else {

      p= ggplot2::ggplot(ggDF,
                         ggplot2::aes(x= Score,
                                      y= Predicted)) +
        ggplot2::geom_point(size= 1.1 + (0.2 * 35/length(model$predicted)),
                            color= "black",
                            alpha= 0.8)
    }


    if(binomial){

      p= p +
        ggplot2::geom_hline(yintercept = 0, color= "dark gray",
                            size= 1.1,
                            linetype= "dashed") +
        ggplot2::ylim(lim) +
        xlab("Group") +
        commonTheme


    } else {

      p= p +
        ggplot2::geom_abline(color= "dark gray",
                             size= 1.1,
                             linetype= "dashed") +
        ggplot2::geom_smooth(method= "lm",
                             color= "red",
                             size= 1.3) +
        ggplot2::ylim(lim) + ggplot2::xlim(lim) +
        commonTheme

    }


  } # end if output== predictions


  if(output== "coefficients"){

    ggDF= model$coeffs

    ord= ggDF$Feature

    ggDF$Feature= as.factor(ggDF$Feature)

    #damn relevel
    for (i in rev(ord)) {

      ggDF$Feature= relevel(ggDF$Feature, i)

    }

    #subset coeffs
    if(!is.null(subset_coeffs)){

      ggDF= ggDF[ggDF$Feature %in% levels(ggDF$Feature)[subset_coeffs],]
      ggDF$Feature= factor(ggDF$Feature)

    }

    #for the fill color
    ggDFt= ggDF$Coefficient

    ggDF$fill_col= abs(ggDFt)

    #plot
    p= ggplot2::ggplot(ggDF,
                       ggplot2::aes(x= Feature,
                                    y= Coefficient,
                                    #group= Feature,
                                    fill= fill_col)) +
      ggplot2::geom_hline(yintercept = 0,
                          color= "dark gray",
                          size= 1.1,
                          linetype= "dashed") +
      ggplot2::scale_fill_gradient(low = "white",
                                   high= "red",
                                   guide= F) +
      commonTheme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         hjust = 1))

    #barplot here
    p= p +
        ggplot2::geom_bar(stat= "identity")


  } #end if output == "coefficients

  #return the plot
  return(p)

}
