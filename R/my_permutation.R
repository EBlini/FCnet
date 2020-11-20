my_permutation= function(y, #dependent variable, typically behavior
                         x, #independent variables, tipycally neural measures
                         alpha= seq(0, 1, by= 0.1),
                         lambda= 10^seq(-5, 5, length.out = 200),
                         nperm= 100, #number of permutations 
                         LOO= F #whether to obtain values from LOO loops or not
){
  
  perms= future_lapply(1:nperm, function(p){
    
    x= my_shuffle(x)
    
    if(LOO){
      
      my_loo(y= y, x= x, alpha = alpha, lambda= lambda)
      
    } else {
      
      my_en(y= y, x= x, alpha = alpha, lambda= lambda)
    }
  })
}
