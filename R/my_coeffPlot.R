my_coeffPlot= function(coeffs){
  
  ggDF= coeffs
  
  ggDFt= abs(tapply(ggDF$Coefficient, 
                    ggDF$Feature, function(x)mean((na.omit(x)))))
  
  #not tvalue anymore
  ggDF$t.value= ggDFt[ggDF$Feature]
  
  #☺this is to show deviant subs on the plot - but it gets cluttered
  
  # g <- ggplot(ggDF, aes(x= Feature, y= Coefficient,
  #                label= Subject,
  #                group= Feature, fill= t.value)) + geom_boxplot()
  # 
  # # get list of outliers 
  # out <- ggplot_build(g)[["data"]][[1]][["outliers"]]
  # 
  # # label list elements with factor levels
  # names(out) <- levels(factor(ggDF$Feature))
  # 
  # tidyout <- map_df(out, as_tibble, .id = "Feature")
  # 
  # foo= function(Feature, value){
  #   
  #   sub= ggDF$Subject[ggDF$Feature == Feature & ggDF$Coefficient== value]
  #   return(sub)
  # }
  # foo= Vectorize(foo)
  # 
  # tidyout$Subject= foo(tidyout$Feature, tidyout$value)
  # 
  # tidyout$Subject= as.character(tidyout$Subject)
  # tidyout$Feature= as.numeric(as.character(tidyout$Feature))
  
  p= ggplot(ggDF, aes(x= Feature, y= Coefficient,
                      label= Subject,
                      group= Feature, fill= t.value)) +
    geom_hline(yintercept = 0, color= "black", size= 1.1, 
               linetype= "dashed") +
    #geom_violin() +
    scale_fill_gradient(low = "white", high= "red",guide= F) +
    geom_boxplot(outlier.colour = "black") +
    # geom_text(data = tidyout, 
    #           aes(x= Feature, y= value,
    #               label = (Subject)), inherit.aes = F,  
    #           hjust = -.3, 
    #           check_overlap = F) +
    ggtitle("Distribution of coefficients")
  
  print(p)
  
  return(p)
}
