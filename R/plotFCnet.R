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
#' @param style For coefficients, how to represent variability (e.g. boxplots or violins).
#' @param plot_labels if `TRUE` annotates the row names of x as predictions.
#' @param subset_coeffs To avoid padding, sometimes the user may want to subset the coefficients to report.
#' This can be done by passing a vector to this parameter. Note that the intercept counts as
#' a coefficient (of 0) even if omitted in the call.
#' @return A `ggplot2` object which can be customized further.

#' @export
#'

plotFCnet= function(model,
                    output= c("predictions", "coefficients"),
                    style= c("boxplot", "violin"),
                    plot_labels= T,
                    subset_coeffs= NULL){


  output= match.arg(output)
  style= match.arg(style)

  #commonTheme
  commonTheme= ggplot2::theme(
    text= ggplot2::element_text(size= 14, face= "bold"),
    axis.text= ggplot2::element_text(size= 14, face= "bold", color= "black"))


  if(output== "predictions"){

    ggDF= data.frame(ID= 1:length(model$predicted),
                     Score= model$y,
                     Predicted= model$predicted)

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


    p= p +
      ggplot2::geom_abline(color= "dark gray",
                           size= 1.1,
                           linetype= "dashed") +
      ggplot2::geom_smooth(method= "lm",
                           color= "red",
                           size= 1.3) +
      commonTheme


  } # end if output== predictions


  if(output== "coefficients"){

    ggDF= model$coeffs

    ord= ggDF$Feature[ggDF$ID== 1]

    ggDF$Feature= as.factor(ggDF$Feature)

    #damn relevel
    for (i in rev(ord)) {

      ggDF$Feature= relevel(ggDF$Feature, i)

    }

    #subset coeffs
    if(!is.null(subset_coeffs)){

      ggDF= ggDF[as.character(ggDF$Feature) %in% subset_coeffs,]
      ggDF$Feature= factor(ggDF$Feature)

    }

    #for the fill color
    ggDFt= abs(tapply(ggDF$Coefficient,
                      ggDF$Feature,
                      function(x)optionsFCnet("consensus_function")((na.omit(x)))))

    ggDF$fill_col= ggDFt[as.character(ggDF$Feature)]


    #plot
    p= ggplot2::ggplot(ggDF,
                       ggplot2::aes(x= Feature,
                                    y= Coefficient,
                                    label= ID,
                                    group= Feature,
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

    if(style== "boxplot"){

      p= p +
        ggplot2::geom_boxplot(outlier.colour = "black")

    } else {

      p= p +
        ggplot2::geom_violin()

    }



  } #end if output == "coefficients

  #return the plot
  return(p)

}