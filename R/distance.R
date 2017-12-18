#' Distance function
#'
#' This function allows you to compute the distance between two curves with the chosen metric.
#' @param FD1 a functional data object of type \code{funData} for the first curve
#' @param FD2 a functional data object of type \code{funData} for the second curve
#' @param metric the chosen distance to be used: \code{"L2"} for the classical L2-distance, \code{"trunc"} for the truncated Mahalanobis semi-distance, \code{"mahalanobis"} for the generalized Mahalanobis distance.
#' @param p a positive numeric value containing the parameter of the regularizing function for the generalized Mahalanobis distance.
#' @param lambda a vector containing the eigenvalues in descending order of the functional data from which the curves are extracted.
#' @param phi a matrix containing the eigenfunctions of the functional data in its columns from which the curves are extracted.
#' @param k_trunc a positive numeric value representing the number of components at which the truncated mahalanobis distance must be truncated
#' @return The function returns a numeric value indicating the distance between the two curves.
#' @keywords distance
#' @references
#'
#' Ghiglietti A., Ieva F., Paganoni A. M. (2017). Statistical inference for stochastic processes:
#' Two-sample hypothesis tests, \emph{Journal of Statistical Planning and Inference}, 180:49-68.
#'
#' Ghiglietti A., Paganoni A. M. (2017). Exact tests for the means of gaussian stochastic processes.
#' \emph{Statistics & Probability Letters}, 131:102--107.
#'
#' @export
#' @examples
#' # Define parameters:
#' n <- 50
#' P <- 100
#' K <- 150
#'
#' # Grid of the functional dataset
#' t <- seq( 0, 1, length.out = P )
#'
#' # Define the means and the parameters to use in the simulation
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
#' z <- gmfd_simulate( n, m1, rho = rho, theta = theta )
#'
#' # Extract two rows of the functional data
#' x <- funData( t, z[1, ] )
#' y <- funData( t, z[2, ] )
#'
#' lambda <- eigen(cov(z))$values
#' phi <- eigen(cov(z))$vectors
#'
#' d <- funDist( x, y, metric = "mahalanobis", p = 1, lambda = lambda, phi = phi )

funDist <- function ( FD1, FD2, metric, p = NULL, lambda = NULL, phi = NULL, k_trunc = NULL ) {
  grid <- FD1$grid
  x <- FD1$data
  y <- FD2$data
  if (!is.list(x)) {
    x <- list(x)
  }
  if (!is.list(y)) {
    y <- list(y)
  }
  R <- length(x) #number of components
  P <- length(grid)
  n <- length(lambda)/R
  # L^2 distance between 2 functions
  if ( metric == "L2" ) {
    Dist <- 0
    for (r in 1:R) {
      integr <- ( x[[r]] - y[[r]] )^2
      Dist <- Dist + integral( grid, integr )
    }
    sqrt(Dist)
  }
  # Truncated Mahalanobis semi-distance
  else if ( metric == "trunc") {
    if ( is.null( lambda ) ) {
      stop( "Eigenvalues are missing!" )
    }
    if ( is.null( phi ) ) {
      stop( "Eigenfunctions are missing!" )
    }
    Dist <- 0
    for ( r in 1:R ) {
      for ( j in 1:k_trunc ) {
        Dist <- Dist + ( integral( grid, ( ( x[[r]] - y[[r]] ) * phi[( 1 + P * ( r - 1 ) ):( r * P ), j] ) ) )^2 / lambda[j]
      }
    }
    Dist <- sqrt( Dist )
  }
  # Generalized Mahalanobis Distance
  else if ( metric == "mahalanobis") {
    h <- NULL
    h <- 1 / (1 / p + lambda)
    if ( is.null( lambda ) ) {
      stop( "Eigenfunctions are missing!" )
    }
    if ( is.null( phi ) ) {
      stop( "Eigenvalues are missing!" )
    }
    Dist <- 0
    for ( r in 1:R ) {
      if ( n < P ) {
        for ( j in 1:( n - 1 ) ) {
          Dist <- Dist + ( integral( grid, ( ( x[[r]] - y[[r]] ) * phi[( 1 + P * ( r - 1 ) ):( r * P ), j] ) ) )^2 * h[j]
        }
        for ( j in n:( P ) ) {
          Dist <- Dist + p * ( integral( grid, ( ( x[[r]] - y[[r]] ) * phi[( 1 + P * ( r - 1 ) ):( r * P ), j] ) ) )^2
        }
      }
      else {
        for (j in 1:( R * P ) ) {
          Dist <- Dist + ( integral( grid, ( ( x[[r]] - y[[r]] ) * phi[( 1 + P * ( r - 1 ) ):( r * P ), j] ) ) )^2 * h[j]
        }
      }
    }
    Dist <- sqrt( Dist )
  }
  return( Dist = Dist )
}


#' Dissimilarity matrix function
#'
#' This function computes the dissimilarity matrix containing the distances between the curves of the functional dataset
#' @param FD a functional data object of type \code{funData}
#' @param metric the chosen distance to be used. Choose \code{"L2"} for the classical L2-distance, \code{"trunc"} for the truncated Mahalanobis semi-distance, \code{"mahalanobis"} for the generalized Mahalanobis distance.
#' @param p a positive numeric value containing the parameter of the regularizing function for the generalized Mahalanobis distance.
#' @param k_trunc a positive numeric value representing the number of components at which the truncated mahalanobis distance must be truncated
#' @return The function returns a matrix of numeric values containing the distances between the curves.
#' @references
#'
#' Ghiglietti A., Ieva F., Paganoni A. M. (2017). Statistical inference for stochastic processes:
#' Two-sample hypothesis tests, \emph{Journal of Statistical Planning and Inference}, 180:49-68.
#'
#' Ghiglietti A., Paganoni A. M. (2017). Exact tests for the means of gaussian stochastic processes.
#' \emph{Statistics & Probability Letters}, 131:102--107.
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
#'
#' FD <- funData( t, x )
#'
#' D <- gmfd_diss( FD, metric = "L2" )

gmfd_diss <- function( FD, metric, p = NULL, k_trunc = NULL ) {
  grid <- FD$grid
  data <- FD$data
  # if data is a univariate functional data, transform it to lists
  if ( !is.list( data ) ) {
    data <- list( data )
  }
  R <- length( data ) #number of components
  if ( R > 1 ) {
    for ( r in 1:( R - 1 ) ) {
      if ( dim( data[[r]] )[[1]] != dim( data[[r + 1]] )[[1]]) {
        stop("All the elements of data must have the same number of rows.")
      }
    }
  }
  P <- length( grid )
  n <- dim( data[[1]] )[[1]]
  Dist <- matrix( 0, n, n )
  h <- NULL
  cluster <- numeric( n )
  xR <- NULL
  for ( r in 1:R )
    xR <- cbind( xR, data[[r]] )
  pca <- eigen( cov( xR ) )
  lambda_hat <- pca$values / ( P * R )
  phi_hat <- pca$vectors * sqrt( ( P * R ) )
  h <- hhat( p, lambda_hat )
  cluster.old <- cluster
  if ( metric == "L2" ) {
    for ( i in 1:n ) {
      for ( j in 1:n ) {
        Dist[i, j] <- sqrt( funDist( funData(grid, lapply( data, '[', i, TRUE )), funData(grid, lapply( data, '[', j, TRUE )), metric = "L2" ) )
      }
    }
  }
  else if ( metric == "trunc" ) {
    for ( i in 1:n ) {
      for ( j in 1:n ) {
        Dist[i, j] <- sqrt( funDist( funData(grid, lapply( data, '[', i, TRUE )), funData(grid, lapply( data, '[', j, TRUE )), metric = "trunc", lambda = lambda_hat, phi = phi_hat, k_trunc = k_trunc ) )
      }
    }
  }
  else if ( metric == "mahalanobis" ) {
    for ( i in 1:n ) {
      for ( j in 1:n ) {
        Dist[i, j] <- sqrt( funDist( funData(grid, lapply( data, '[', i, TRUE )), funData(grid, lapply( data, '[', j, TRUE )), metric = "mahalanobis", lambda = lambda_hat, phi = phi_hat, p = p ) )
      }
    }
  }
  return( Dist = Dist )
}
