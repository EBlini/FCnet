my_shuffle= function(x){
  
  x= x[sample(1:nrow(x), size = nrow(x), replace = F),]
  
  return(x)
  
}