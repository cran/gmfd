#' k-means clustering algorithm
#'
#' This function performs a k-means clustering algorithm on an univariate or multivariate functional data using a generalization of Mahalanobis distance.
#' @param FD a functional data object of type \code{funData}.
#' @param n.cl an integer representing the number of clusters.
#' @param metric the chosen distance to be used: \code{"L2"} for the classical L2-distance, \code{"trunc"} for the truncated Mahalanobis semi-distance, \code{"mahalanobis"} for the generalized Mahalanobis distance.
#' @param p a positive numeric value containing the parameter of the regularizing function for the generalized Mahalanobis distance.
#' @param k_trunc a positive numeric value representing the number of components at which the truncated mahalanobis distance must be truncated
#' @keywords Clustering
#' @return The function returns a list with the following components:
#' \code{cluster}: a vector of integers (from \code{1} to \code{n.cl}) indicating the cluster to which each curve is allocated;
#' \code{centers}: a list of \code{d} matrices (\code{k} x \code{T}) containing the centroids of the clusters
#' @references
#'
#' Martino A., Ghiglietti A., Ieva F., Paganoni A. M. (2017). A k-means procedure based on a Mahalanobis type
#' distance for clustering multivariate functional data, \emph{MOX report 44/2017}
#'
#' Ghiglietti A., Ieva F., Paganoni A. M. (2017). Statistical inference for stochastic processes:
#' Two-sample hypothesis tests, \emph{Journal of Statistical Planning and Inference}, 180:49-68.
#'
#' Ghiglietti A., Paganoni A. M. (2017). Exact tests for the means of gaussian stochastic processes.
#' \emph{Statistics & Probability Letters}, 131:102--107.
#'
#' @seealso \code{\link{funDist}}
#' @import
#' graphics
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
#' m1 <- t^2 * ( 1 - t )
#'
#' rho <- rep( 0, K )
#' theta <- matrix( 0, K, P )
#' for ( k in 1:K) {
#'   rho[k] <- 1 / ( k + 1 )^2
#'   if ( k%%2 == 0 )
#'     theta[k, ] <- sqrt( 2 ) * sin( k * pi * t )
#'   else if ( k%%2 != 0 && k != 1 )
#'     theta[k, ] <- sqrt( 2 ) * cos( ( k - 1 ) * pi * t )
#'   else
#'     theta[k, ] <- rep( 1, P )
#' }
#'
#' s <- 0
#' for (k in 4:K) {
#'  s <- s + sqrt( rho[k] ) * theta[k, ]
#' }
#'
#' m2 <- m1 + s
#'
#' # Simulate the functional data
#' x1 <- gmfd_simulate( n, m1, rho = rho, theta = theta )
#' x2 <- gmfd_simulate( n, m2, rho = rho, theta = theta )
#'
#' # Create a single functional dataset containing the simulated datasets:
#' FD <- funData(t, rbind( x1, x2 ) )
#'
#' output <- gmfd_kmeans( FD, n.cl = 2, metric = "mahalanobis", p = 10^6 )

gmfd_kmeans <- function( FD, n.cl = 2, metric, p = NULL, k_trunc = NULL ) {
  grid <- FD$grid
  data <- FD$data
  if ( n.cl - round( n.cl ) != 0 ) {
    stop( "n.cl must be an integer number." )
  }
  R <- length( data ) #number of components
  if ( R > 1 ) {
    for ( r in 1:( R - 1 ) ) {
      if ( dim( data[[r]] )[[1]] != dim( data[[r + 1]] )[[1]] ) {
        stop( "All the elements of data must have the same number of rows." )
      }
    }
  }
  P <- length( grid )
  n <- dim( data[[1]] )[[1]]
  C <- list( )
  s <- sample( n, n.cl )
  for ( r in 1:R ) {
    C[[r]] <- data[[r]][s, ] #centroids of indices s of the r-th functions
  }
  Dist <- matrix( 0, n, n.cl )
  cluster <- rep( 1, n )
  cluster.old <- rep( 2, n )
  cluster <- numeric(n)
  xR <- NULL
  for ( r in 1:R )
    xR <- cbind( xR, data[[r]] )
  pca <- eigen( cov( xR ) )
  lambda_hat <- pca$values / ( P*R )
  phi_hat <- pca$vectors * sqrt( ( P*R ) )
  while( sum( abs( cluster - cluster.old ) ) > 0) {
    cluster.old <- cluster
    if (metric == "L2") {
      for (i in 1:n) {
        for (k in 1:n.cl){
          Dist[i, k] <- ( funDist( funData(grid, lapply( data, '[', i, TRUE )), funData(grid, lapply( C, '[', k, TRUE )), metric = "L2" ) )
        }
      }
    }
    else if ( metric == "trunc" ) {
      for ( i in 1:n ) {
        for ( k in 1:n.cl ) {
          Dist[i, k] <- ( funDist( funData(grid, lapply( data, '[', i, TRUE )), funData(grid, lapply( C, '[', k, TRUE )), metric = "trunc", lambda = lambda_hat, phi = phi_hat, k_trunc = k_trunc ) )
        }
      }
    }

    else if ( metric == "mahalanobis" ) {
      for ( i in 1:n ) {
        for ( k in 1:n.cl ) {
          Dist[i, k] <- ( funDist( funData(grid, lapply( data, '[', i, TRUE )), funData(grid, lapply( C, '[', k, TRUE )), metric = "mahalanobis", lambda = lambda_hat, phi = phi_hat, p = p) )
        }
      }
    }

    for( i in 1:n )
      cluster[i] <- which.min( Dist[i, ] )

    C <- list( )
    for ( r in 1:R ) {
      CC <- NULL
      for( k in 1:n.cl ) {
        CC <- rbind( CC, colMeans( data[[r]][cluster == k, ] ) )
      }
      C[[r]] <- CC
    }
  }
  par( mfrow = c(1, 2), ask = F)
  for ( r in 1:R ) {
    plot( grid, data[[r]][1, ], type = "l", col = 2, lty = 3, main = paste('Groups for component ', r), ylim = c( min( data[[r]] ) - 2, max( data[[r]]) + 2 ), xlab = "t", ylab = "X(t)" )
    for ( i in 1:n ) {
      for ( k in 1:n.cl ) {
        if( cluster[i] == k )
          lines( grid, data[[r]][i, ], col = 1 + k, lty = 3 )
      }
    }
    matplot( grid, t( data[[r]] ), type = "l", col = "grey", main = paste( 'Centers for component ', r ), lty = 3, xlab = "t", ylab="X(t)" )
    for ( k in 1:n.cl ) {
      lines( grid, C[[r]][k, ], col = 1 + k, lwd = 3 )
    }
    par( mfrow = c( 1, 2 ) )
  }
  list( cluster = cluster,
        centers = C )
}
