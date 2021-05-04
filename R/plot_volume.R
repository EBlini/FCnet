#' plot a brain Volume
#'
#' This function reads a 3-D array and returns a very basic plot
#' using `neurobase::ortho2()`. This can be, for example, a back-projected volume
#' or the mean volume returned from the reduce_features object.
#'
#' @param volume A 3_d array to depict.
#' @param x Coordinates at which the orthostatic view is displayed.
#' @param y Coordinates at which the orthostatic view is displayed.
#' @param z Coordinates at which the orthostatic view is displayed.
#' @param limit Range of values to depict. If NULL, the default,
#' simmetry is preserved around the largest value in absolute value.
#' @param col.y The color palette used to plot values (i.e. a vector of colors). The vector is often divergent, es. `c("red", "white", "blue")`.
#' Defaults to `optionsFCnet("colorPaletteDefault")`.
#' @param colorNuances Number of nuances along the provided col.y.

#' @export

plot_volume= function(volume,
                      x= 1,
                      y= 1,
                      z= 1,
                      limit= NULL,
                      col.y= c("darkslateblue", "white", "firebrick4"),
                      colorNuances= optionsFCnet("colorNuances")){

  #"smart" limit setting
  if(is.null(limit)){

    min.l= min(volume, na.rm = T)
    max.l= max(volume, na.rm = T)

    if(abs(min.l)<abs(max.l)){

      min.l= -1*max.l

    } else {

      max.l= -1*min.l

    }

    limit= c(min.l, max.l)

  }

  #empty array
  zeros= base::array(0, dim(volume))

  #get colors
  get_col= colorRampPalette(col.y)
  col_y= get_col(colorNuances-1)

  #print(
  neurobase::ortho2(zeros, volume,
         xyz = c(x, y, z),
         col.y = col_y,
         crosshairs = F,
         addlegend = T,
         legend= "",
         leg.cex = 1.5,
         leg.col = "black",
         ybreaks = seq(limit[1], limit[2],
                       length.out = colorNuances),
         ycolorbar = T)
  #)



}







