#' Obtain the vectorisation of a matrix
#'
#' Given an \code{n x m} matrix, returns an \code{(nm)\times 1} matrix
#' containing the column
#' stack of the matrix
#' @param A matrix to vectorize
#' @return one-column matrix; the vectorization of A
#' @seealso \code{\link{vecInverse}}, \code{\link{as.vector}}
#' @export
vec <- function(A){
    return(as.vector(A))
}
