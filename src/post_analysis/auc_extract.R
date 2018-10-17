library(ROCR)

# it extracts off-diagonal elements of an matrix
# returns a vector
offdiag.extract <- function(m) {
  m_upper = m[upper.tri(m)]
  m_lower = m[lower.tri(m)]
  
  return(c(m_upper, m_lower))
}

block.extract <- function(m, block_bool){
  return(m[block_bool])
}

auc.extract <- function(y.pred, y.lab){
  # return the numerical value
  y.pred <- offdiag.extract(y.pred)
  y.lab <- offdiag.extract(y.lab)
  
  pred <- prediction(as.numeric(y.pred), as.numeric(y.lab))
  return(performance(pred, "auc")@y.values[[1]])
}
