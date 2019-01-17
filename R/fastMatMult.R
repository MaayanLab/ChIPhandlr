#' @export

myFastMult = function(A,B){
  return(eigenMapMatMult(A,B))
}

#' @export
multByCols = function(A,B){
  for(i in 1:length(B)){
    A[,i] = A[,i]*B[i]
  }
  return(A)
}