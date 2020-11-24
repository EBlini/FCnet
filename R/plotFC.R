#' plot a Functional Connectivity matrix
#'
#' This function reads a (square) matrix and plots it using `ggplot2`.
#'
#' @param FCmatrix The FC matrix to plot
#' @param style Whether to plot only the lower triangle (and diagonal) or the entire FC matrix.
#' @param col The color palette used to plot values (i.e. a vector of colors). The vector is often divergent, es. `c("red", "white", "blue")`.
#' Defaults to `optionsFCnet("colorPaletteDefault")`.
#' @param limit The limits for the FC values. Defaults to `NULL` and automatically adapts to the data range.
#' @param network_definition A character vector specifying to which FC network the ROI belongs to. If provided, draws vertical and horizontal lines visually separating the networks.
#' @param plot_labels if `TRUE` annotates the network names defined in `network_definition`. Pretty results not warranted.
#'
#' @return A `ggplot2` object which can be customized further.

#' @export

plotFC= function(FCmatrix,
                 style= c("lower.tri", "full"),
                 col= optionsFCnet("colorPaletteDefault"),
                 limit = NULL,
                 network_definition= NULL,
                 plot_labels= F){

  #match
  style= match.arg(style)

  if(style== "lower.tri"){

    #erase upper triangle
    FCmatrix[lower.tri(FCmatrix,
                       diag = F)]= NA

    }

  #matrix to long dataframe
  hmDF= reshape2::melt(as.matrix(FCmatrix))

  # #if inf are there...
  hmDF$value[hmDF$value== Inf]= NA
  hmDF$value[hmDF$value== -Inf]= NA

  if(is.null(limit))(limit= c(min(hmDF$value, na.rm = T),
                              max(hmDF$value, na.rm = T)))

  #dimension of the matrix
  mat_dim= dim(FCmatrix)[1]

  #rename levels
  levels(hmDF[,2])= 1:mat_dim
  hmDF[,2]= as.numeric(hmDF[,2])

  #retrieve color scale from options
  gc= colorRampPalette(col)
  col= gc(20)

  #plot baseline grid
  p= ggplot2::ggplot(data = hmDF,
                     ggplot2::aes(x=Var1, y=Var2, fill=value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = col,
                         limit = limit,
                         na.value= "transparent") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::scale_y_continuous(trans = "reverse") +
    ggplot2::guides(fill= NULL) +
    ggplot2::theme(axis.title.x= ggplot2::element_blank(),
          axis.text.x= ggplot2::element_blank(),
          axis.title.y= ggplot2::element_blank(),
          axis.text.y= ggplot2::element_blank())


  #add network labels if provided
  if(!is.null(network_definition)){

    #intercepts
    int= rle(network_definition)

    intercepts= data.frame(intercepts= c(0, cumsum(int$lengths)),
                           labels= c(NA, int$values))

    #Dintercepts coordinates change as a function of style
    if(style== "lower.tri"){
      int_v= intercepts$intercepts
      int_v2= rep(mat_dim,
                  length(intercepts$intercepts))
      } else {
        int_v= rep(mat_dim,
                   length(intercepts$intercepts))
        int_v2= rep(0,
                    length(intercepts$intercepts))
      }

    p=  p + ggplot2::geom_segment(
        data= data.frame(intercepts= intercepts$intercepts,
                         x= rep(0,
                                length(intercepts$intercepts)),
                         xend= int_v
                         ),
        ggplot2::aes(y= intercepts, yend= intercepts, x= x, xend= xend),
        inherit.aes = F) +
      ggplot2::geom_segment(
        data= data.frame(intercepts= intercepts$intercepts,
                         yend= int_v2,
                         y= int_v
                         ),
        ggplot2::aes(y= y, yend= yend, x= intercepts, xend= intercepts),
        inherit.aes = F)

    #plot diagonal anyway
    p= p + ggplot2::geom_segment(
      data= data.frame(x= 0,
                       y= 0,
                       xend= mat_dim,
                       yend= mat_dim),
      ggplot2::aes(y= y, yend= yend, x= x, xend= xend), inherit.aes = F)



    if(plot_labels){

      labels= data.frame(
        pos= rowMeans(embed(intercepts$intercepts, 2)),
        labels= int$values)

      if(style== "lower.tri"){
        #check spacing
        #you may bring this function out
        cs= function(x){

          while(sum(diff(x) < mean(diff(x)) - 1.6*sd(diff(x)))>0){
            co= mean(diff(x)) - 1.6*sd(diff(x))
            dpos= c(FALSE, diff(x)< co)
            #x[which(dpos)-1]= x[which(dpos)-1] - 0.1
            x[which(dpos)]= x[which(dpos)] + 0.1

          }

          return(x)
        }
        labels$pos= cs(labels$pos)

        p= p +
          ggplot2::annotate("text", x = labels$pos+1,
                   y= labels$pos,
                   label= labels$labels,
                   hjust= 0, vjust= 0,
                   fontface= "bold")
      } else {

      labels$alt= rep(c(-6, mat_dim+5,
                        -14, mat_dim+13),
                      200)[1:length(labels$pos)]

      labels$angle= rep(c(0, 0),
                        200)[1:length(labels$pos)]

      labels$alt2= rep(c(-6, mat_dim+5,
                         -14, mat_dim+13),
                       200)[1:length(labels$pos)]

      labels$angle2= rep(c(90, -90),
                         200)[1:length(labels$pos)]

      p= p +
        ggplot2::annotate("text", x = labels$pos,
                 y= labels$alt,
                 label= labels$labels,
                 angle= labels$angle, hjust= 0.5) +
        ggplot2::annotate("text", y = labels$pos,
                 x= labels$alt2,
                 label= labels$labels,
                 angle= labels$angle2, hjust= 0.5)

      } #end else switch

    } #end if plot_labels

  } #end if network is defined

  #print(p)
  return(p)
}
