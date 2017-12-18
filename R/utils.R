integral <- function ( grid, values ) {
  return( sum( ( values[ -1 ] + values[ -length( values ) ] ) * diff( grid ) / 2 ) )
}


hhat <- function( p , lambda ) {  	# definisco la funzione h
  return( lambda / ( lambda + 1 / p ) )
}
