ramsvm <- function(x, 
                   y, 
                   lambda, 
                   gamma = 0.5, 
                   weight = NULL, 
                   kernel = "linear", 
                   kparam = NULL, 
                   large = FALSE, 
                   epsilon = NULL, 
                   warm = NULL, 
                   nb.core = NULL) {

  #------------------------------------------------------------------#
  # Verify that kernel is one of linear, polynomial, gaussian        #
  #------------------------------------------------------------------#
  kernel <- tolower(x = kernel)
  if( !(kernel %in% c("linear","polynomial", "gaussian")) ) {
    stop("kernel must be one of ('linear', 'polynomial', 'gaussian')")
  }

  if( kernel %in% c("polynomial", "gaussian") ) {
    if( is.null(x = kparam) ) kparam = 1.0
  }

  #------------------------------------------------------------------#
  # Verify that covariates are provided in matrix form with no NA/NAN#
  #------------------------------------------------------------------#
  if( !is.matrix(x = x) ) stop("The covariates must be a matrix.")  
  if( any(is.na(x = x)) || any(is.nan(x = x)) ) {
    stop("There must be no NA/NaN in the covariates.")
  }

  nobs <- nrow(x)

  #------------------------------------------------------------------#
  # Verify that the number of labels matches the number of obs.      #
  #------------------------------------------------------------------#
  if( length(x = y) != nobs ) {
    stop("the dimension of covariates must match the length of the label")
  }

  #------------------------------------------------------------------#
  # determine the number of classes; must be > 1.                    #
  #------------------------------------------------------------------#
  k <- length(x = levels(x = as.factor(x = y)))
  if( k < 2L) stop("there must be at least two classes")
  
  if( is.null(x = epsilon) ) {
    #--------------------------------------------------------------#
    # If not provided, calculate epsilon.                          #
    #--------------------------------------------------------------#
    epsilon <- 0.0001 * nobs * k
  } else {
    #--------------------------------------------------------------#
    # If provided, verify that epsilon is a positive scalar.       #
    #--------------------------------------------------------------#
    if( is.na(x = epsilon) || is.nan(x = epsilon) ) {
      stop("epsilon must not be NA/NaN")
    }
    if( !is.numeric(x = epsilon) ) stop("epsilon must be numeric")
    if( epsilon <= 0.0 ) stop("epsilon must be strictly positive")
    if( length(x = epsilon) > 1.5 ) stop("epsilon must be a scalar")
  }
    
  if( is.null(x = weight) ) {
    #--------------------------------------------------------------#
    # If not provided, set weight to 1.0                           #
    #--------------------------------------------------------------#
    weight <- numeric(nobs) + 1.0
  } else {
    #--------------------------------------------------------------#
    # If provided, verify that weight is a positive, numeric vector#
    # of length equivalent to the number of observations.          #
    #--------------------------------------------------------------#
    if( !is.numeric(x = weight) ) {
      stop("weight vector must be numeric")
    }
    if( any(is.na(x = weight)) || any(is.nan(x = weight)) ) {
      stop("no NA/NaN allowed in weight vector")
    }
    if( min(weight) < 0.0 ) {
      stop("weight vector must be non-negative")
    }
    if( length(x = weight) != nobs ) {
      stop("length of weight must be the number of observations")
    }
  }

  if( !is.null(x = warm) ) {
    #--------------------------------------------------------------#
    # If provided, verify dims of starter matrix are appropriate.  #
    #--------------------------------------------------------------#
    if( {nrow(x = warm) != nobs} || {ncol(x = warm) != k} ) {
      stop("dimension of the warm start matrix is incorrect")
    }
  } else {
    #--------------------------------------------------------------#
    # If not provided, default starter matrix to zero              #
    #--------------------------------------------------------------#
    warm <- matrix(data = 0.0, nrow = nobs, ncol = k)
  }

  #------------------------------------------------------------------#
  # Verify that all lambdas are positive.                            #
  #------------------------------------------------------------------#
  if( any(lambda <= 0) ) stop("all lambdas must be positive.")

  #------------------------------------------------------------------#
  # If the lambdas are not sorted warn user that they are reorder.   #
  #------------------------------------------------------------------#
  if( any( order(lambda) != order(sort(lambda,TRUE))) ) {
    warning("order of lambda has been changed")
  }

  #------------------------------------------------------------------#
  # Sort lambda high to low.                                         #
  #------------------------------------------------------------------#
  lambda = as.double(rev(sort(lambda)))

  if( is.null(x = nb.core) && large ) {
    #--------------------------------------------------------------#
    # If this is a large calculation and the number of cores are   #
    # not specified, detect the number of available cores.         #
    #--------------------------------------------------------------#
    nb.core <- parallel::detectCores()
  } else if( !large ) {
    nb.core <- 0L
  }

  #------------------------------------------------------------------#
  # Call solver method.                                              #
  #------------------------------------------------------------------#
  fit <- RAMSVM_solve(x = x,
                      y = y,
                      gamma = gamma,
                      lambda = lambda,
                      kernel = kernel,
                      kparam = kparam,
                      weight = weight,
                      epsilon = epsilon,
                      warm = warm,
                      nb.core = nb.core)

  fit@call <- match.call()
  return(fit)
}
