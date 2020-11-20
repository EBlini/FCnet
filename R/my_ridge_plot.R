my_ridge_plot= function(Subject, score, predicted, name){
  
  ggDF= data.frame(Subject= Subject,
                   Score= score,
                   Predicted= predicted)
  
  p= ggplot(ggDF, aes(x= Score, y= Predicted, label= Subject)) +
    geom_label() + 
    geom_abline(color= "black", size= 1.1, linetype= "dashed") + 
    geom_smooth(method= "lm", color= "red", size= 1.3) +
    ggtitle(name) 
  
  #p= recordPlot()
  print(p)
  return(p)
  
}

