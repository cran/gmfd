#' Simulation of a functional sample
#'
#' Simulate a univariate functional sample using a Karhunen Loeve expansion.
#' @param size a positive integer indicating the size of the functional sample to simulate.
#' @param mean a vector representing the mean of the sample.
#' @param covariance a matrix from which the eigenvalues and eigenfunctions must be extracted.
#' @param rho a vector of the eigenvalues in descending order to be used for the simulation.
#' @param theta a matrix containing the eigenfunctions in its columns to be used for the simulation.
#' @keywords Simulation
#' @return The function returns a functional data object of type \code{funData}.
#'
#' @import
#' stats
#' @export
#' @examples
#' # Define parameters
#' n <- 50
#' P <- 100
#' K <- 150
#'
#' # Grid of the functional dataset
#' t <- seq( 0, 1, length.out = P )
#'
#' # Define the means and the parameters to use in the simulation
#' # with the Karhunen - Loève expansion
#' m1 <- t^2 * ( 1 - t )
#'
#' rho <- rep( 0, K )
#' theta <- matrix( 0, K, P )
#' for ( k in 1:K ) {
#'   rho[k] <- 1 / ( k + 1 )^2
#'   if ( k%%2 == 0 )
#'     theta[k, ] <- sqrt( 2 ) * sin( k * pi * t )
#'   else if ( k%%2 != 0 && k != 1 )
#'     theta[k, ] <- sqrt( 2 ) * cos( ( k - 1 ) * pi * t )
#'   else
#'     theta[k, ] <- rep( 1, P )
#' }
#'
#' # Simulate the functional data
#' x <- gmfd_simulate( n, m1, rho = rho, theta = theta )


gmfd_simulate <- function( size, mean, covariance = NULL, rho = NULL, theta = NULL ) {
  P <- length( mean )
  u <- matrix( 0, size, P )
  if( is.null( covariance ) ) {
    K_trunc <- length( rho )
  }
  else {
    K_trunc <- dim( covariance )[2]
  }
  z <- matrix( 0, size, K_trunc )
  if ( is.null( rho ) | is.null( theta ) ) {
    eig <- eigen( covariance )
    rho <- eig$values
    theta <- eig$vectors
  }
  for ( i in 1:size ) {
    z[i, ] <- rnorm( K_trunc )
    for( k in 1:K_trunc ) {
      u[i, ] <- u[i, ] + sqrt( rho[k] ) * ( z[i, k] * theta[k, ] )
    }
  }
  x <- matrix( 0, size, P )
  for ( i in 1:size )
    x[i, ] <- u[i, ] + mean
  return( x )
}
