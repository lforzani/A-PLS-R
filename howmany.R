#####################
##function that counts how many times a quantity 
##belongs to a sequence of intervals needed for Table 3  
#####################

howmany <- function(X,y){
  s <- 0
  values <- NULL
  ind <- NULL
  vec <- NULL
  for (k in 1:length(y)){
    if (X[k,1] <= y[k] & y[k]<=X[k,2]){
      s <- s+1
      values[k] <- y[k]
      ind[k] <- k
      vec[k] <- 1
      }
    else{
      s <- s
      vec[k] <- 0
    }  
  }
  ind <- ind[!is.na(ind)]
  values <- values[!is.na(values)]
  percen <- s/length(y)
  out <- list(s,percen,values,ind,vec)
  names(out) <- c("amount","percentage","values","index","vector")
  return(out)
}
