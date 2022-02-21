tobitnet = function(x, y, c = 0, nlambda = 100, lambda.factor = ifelse(n/2 < nvars, 0.05, 0.01), 
                    lambda1 = NULL, lambda2 = 0, pf1 = rep(1, nvars), pf2 = rep(1, nvars), eps = 1e-7, 
                    standardize = TRUE, maxit = 1e6, early.stop = TRUE){
    this.call = match.call()
    
    n = nrow(x)
    nvars = ncol(x)
    
    #x is numeric, a matrix, and has finite entries
    if( any( c(is.null(nvars), nvars <= 1) ) ) stop("x must be a matrix with 2 or more columns")
    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.matrix(x)) stop("x must be a matrix")
    if(!is.numeric(x) | !all(is.finite(x))) stop("x must be a design matrix with only finite, numeric entries")
    
    #Check that no column of x has colsd = 0
    x_check = scale(x)
    constant_cols <- attributes(x_check)$'scaled:scale' == 0
    if( any(constant_cols) ) {
        #Remove constant columns from x
        x = x[, !constant_cols]
        x_check = x_check[, !constant_cols]
        
        #Remove penalty factors corresponding to removed columns
        pf1 = pf1[!constant_cols]
        pf2 = pf2[!constant_cols]
    }
    p = ncol(x)
    
    #y is numeric, a vector, and has finite entries and nrow(x) = length(y)
    if(!is.vector(y)) y <- as.vector(y)
    if(!is.vector(y)) stop("y must be a vector")
    if(length(y) != n) stop("the number of observations in y is not equal to the number of rows in x")
    if(!is.numeric(y) | !all(is.finite(y))) stop("y must be a vector with only finite, numeric entries")
    
    #c is numeric and scalar
    if(!is.numeric(c) | length(c) != 1 | !is.finite(c)) stop("c must be a finite scalar") 
    
    #lambda1 and lambda2 are numeric, positive, and finite
    if(!is.null(lambda1)){
        if(!is.vector(lambda1)) lambda1 <- as.vector(lambda1)
        if(!is.vector(lambda1)) stop("lambda1 must be a vector if it is not set to NULL")
        if(!is.numeric(lambda1) | !all(lambda1 >= 0) | !all( is.finite(lambda1) )) stop("lambda1 must be a vector with only non-negative, finite numeric entries if it is not set to NULL")
    }
    if(!is.numeric(lambda2) | !(lambda2 >= 0) | !( is.finite(lambda2) ) | length(lambda2) != 1 ) stop("lambda2 must be a single non-negative, finite number")
    
    #lambda.factor is numeric, scalar, and < 1
    if(!is.numeric(lambda.factor) | length(lambda.factor) != 1 | !(lambda.factor < 1) | !(lambda.factor > 0) ) stop("lambda.factor must be a scalar between 0 and 1")
    
    #If lambda1 = NULL, nlambda > 2. nlambda must be a positive integer
    if(is.null(lambda1)){
        if(!is.numeric(nlambda) | nlambda != round(nlambda) | !(nlambda > 2) | !( is.finite(nlambda) ) | length(nlambda) != 1 ) stop("nlambda must be a single integer greater than 2")
    }
    
    #pf1 and pf2 are numeric vectors. length(pf1) == length(pf2) == p
    if(!is.vector(pf1)) pf1 <- as.vector(pf1)
    if(!is.vector(pf1)) stop("pf1 must be a vector")
    if(!is.vector(pf2)) pf2 <- as.vector(pf2)
    if(!is.vector(pf2)) stop("pf2 must be a vector")
    stopifnot( length(pf1) == p, length(pf2) == p, is.numeric(pf1), is.numeric(pf2))
    
    #eps is numeric, positive, and finite
    if(!is.numeric(eps) | !(eps > 0) | !( is.finite(eps) ) | length(eps) != 1 ) stop("eps must be a single positive, finite number")
    
    #standardize is boolean
    stopifnot( is.logical(standardize), length(standardize) == 1 )
    
    #maxit must be a positive integer
    if(!is.numeric(maxit) | maxit != round(maxit) | !(maxit > 0) | !( is.finite(maxit) ) | length(maxit) != 1 ) stop("maxit must be a single positive integer")
    
    #early.stop is boolean
    stopifnot( is.logical(early.stop), length(early.stop) == 1 )
    
    d = (y > c)
    
    #Fit null model and get fitted values
    mod_null = tobitnet_innerC(xin = matrix(0, 0, 0), yin = y, cin = c, lambda1 = 0, lambda2 = lambda2, pf1 = pf1, pf2 = pf2, delta_init = rep(0,p), eps = eps, standardize = F, maxit = maxit)
    if(!mod_null$KKT) warning("KKT conditions not satisfied. Try decreasing eps or increasing maxit.")
    r_null = rep(mod_null$d0, n)
    gamma_null = mod_null$gamma
    
    #Set lambda1 solution path if not provided
    if(is.null(lambda1)){
        if(!standardize) x_check = x
        
        #Compute lambda max
        l1max = max( vapply(1:p, function(j) abs( sum( LprimeC(xj = x_check[,j], y = y-c, d = d, r = r_null, gamma = gamma_null) )/n ), FUN.VALUE = numeric(1) ) )
        l1min = l1max*lambda.factor
        
        llseq = seq(log(l1max), log(l1min), length.out = nlambda)
        lambda1 = exp(llseq)
        
        delta_0_init = mod_null$d0
        delta_init = rep(0,p)
        gamma_init = mod_null$gamma
    } else {
        delta_0_init = 0 
        delta_init = rep(0,p)
        gamma_init = 1
        nlambda = length(lambda1)
    }
    
    gamma_vec = rep(0, nlambda)
    b0_vec = rep(0, nlambda)
    beta_mat = matrix(0, nrow = p, ncol = nlambda)
    
    #Compute null deviance
    null_dev = -2*sum( logL1(y = y-c, d = d, r = r_null, gamma = gamma_null) )
    dev = rep(0, nlambda)
    
    #Fit sequence of tobitnets along lambda1 path
    for(l in 1:nlambda){
        l1 = lambda1[l]
        tn = tobitnet_innerC(xin = x, yin = y, cin = c, lambda1 = l1, lambda2 = lambda2, pf1 = pf1, pf2 = pf2,
                            delta_0_init = delta_0_init, delta_init = delta_init, gamma_init = gamma_init,
                            eps = eps, standardize = standardize, maxit = maxit)
        if(!tn$KKT) warning("KKT conditions not satisfied. Try decreasing eps or increasing maxit.")
        gamma_vec[l] = tn$gamma
        b0_vec[l] = tn$b0
        beta_mat[,l] = tn$beta
        
        #Compute deviance
        delta_pred = matrix(tn$beta*tn$gamma, nrow = p, ncol = 1)
        r_temp = x%*%delta_pred + rep((tn$b0 - c)*tn$gamma, n)
        dev[l] = -2*sum( logL1(y = y-c, d = d, r = r_temp, gamma = tn$gamma) )

        #Stop early if the deviance is barely changing
        if(early.stop == T & l > 2){
            if( abs( (dev[l] - dev[l-1])/null_dev ) <= 1e-5 ){
                nlambda = l
                break
            }
        }
        
        #Set previous solution as new inits (Note: tn$d0 and tn$delta are in the scale of the standardized predictors if standardize == T)
        delta_0_init = tn$d0
        delta_init = tn$delta
        gamma_init = tn$gamma
    }
    
    #Put in rows of 0s for betas corresponding to the constant columns
    beta_final = matrix(0, nrow = nvars, ncol = nlambda)
    beta_final[!constant_cols, ] = beta_mat[, 1:nlambda, drop = F]
    
    return(structure( list(
                call = this.call,
                sigma = 1/gamma_vec[1:nlambda],
                b0 = b0_vec[1:nlambda],
                beta = beta_final,
                c = c,
                lambda1 = lambda1[1:nlambda], 
                lambda2 = lambda2,
                dev = dev[1:nlambda],
                nulldev = null_dev
                ), class = "tobitnet"))
}

predict.tobitnet = function(object, newx, lambda1 = NULL, type = c("censored", "uncensored"), ...){
    type = match.arg(type)
    this.call = match.call()
    
    if(missing(newx)) stop("Missing a value for newx")
    #newx is numeric, a matrix, and has finite entries
    if(is.vector(newx)){
      newx = matrix(newx, nrow = 1, ncol = length(newx))
    } else if(!is.matrix(newx)){ 
      newx = as.matrix(newx)
    }
    if(!is.matrix(newx)) stop("newx must be a matrix")
    if(!is.numeric(newx) | !all(is.finite(newx))) stop("newx must be a design matrix with only finite, numeric entries")

    n = nrow(newx)
    p = ncol(newx)
    
    if(p != nrow(object$beta)) stop("newx must have the same number of columns as the matrix x used to fit the model")
    
    c = object$c
    
    if(!is.null(lambda1)){ 
        object = update(object, lambda1 = lambda1, c = c , ...)
    }

    beta_0 = matrix(object$b0, nrow = 1, ncol = length(object$lambda1))
    
    r = newx%*%object$beta + matrix(1, nrow = n, ncol = 1)%*%beta_0
    if(type == 'censored'){
        r = pmax(r, c)
    }
    
    return( r )
}

cv.tobitnet = function(x, y, c = 0, lambda1 = NULL, nfolds = 10, early.stop = TRUE, ...){
    this.call = match.call()
    n = nrow(x)
    p = ncol(x)
    
    if( any( c(is.null(p), p <= 1) ) ) stop("x must be a matrix with 2 or more columns")
    
    #nfolds must be an integer between 2 and n (other argument checks covered by tobitnet)
    if(!is.numeric(nfolds) | nfolds != round(nfolds) | !(nfolds >= 2) | !(nfolds <= n) | !( is.finite(nfolds) ) | length(nfolds) != 1 ) stop("nfolds must be a single positive integer between 2 and the number of observations")
    
    tn_init = tobitnet(x = x, y = y, c = c, lambda1 = lambda1, early.stop = early.stop, ...)
    lambda1 = tn_init$lambda1
    nlambda = length(lambda1)
    
    #Fold assignments
    nonzero_indices = which(y > c)
    n_nz = length(nonzero_indices)

    zero_indices = which(y == c)
    n_z = length(zero_indices)

    nperfold_nz = floor(n_nz/nfolds)
    unfolded_entries_nz = nonzero_indices

    nperfold_z = floor(n_z/nfolds)
    unfolded_entries_z = zero_indices
    
    #CV loop setup
    foldlist <- list()
    err_mat = matrix(0, nrow = nfolds, ncol = nlambda)
    
    for(i in 1:nfolds){
        if(i < nfolds){
            fold_nz = sample(unfolded_entries_nz, nperfold_nz)
            fold_z = sample(unfolded_entries_z, nperfold_z)

            unfolded_entries_nz = setdiff(unfolded_entries_nz, fold_nz)
            unfolded_entries_z = setdiff(unfolded_entries_z, fold_z)
        } else {
            fold_nz = unfolded_entries_nz
            fold_z = unfolded_entries_z
        }
        foldlist[[i]] = c(fold_nz, fold_z)
        
        tn = tobitnet(x = x[-foldlist[[i]], ], y = y[-foldlist[[i]] ], c = c, lambda1 = lambda1, early.stop = F, ...)
        
        #Predict and test
        preds = predict(tn, newx = x[foldlist[[i]],], type = "censored") #n x nlambda
        
        r2 = (matrix(y[ foldlist[[i]] ], nrow = length(foldlist[[i]]), ncol = nlambda) - preds)^2
        err_mat[i,] = colMeans(r2) 
    }
    
    cvm = colMeans(err_mat) 
    cvvar = colMeans(err_mat*err_mat) - cvm^2 
    cvsd = sqrt(cvvar/nfolds)
    lambda1.min = lambda1[which.min(cvm)]
    lambda1.1se = max(lambda1[ cvm < min(cvm) + cvsd[which.min(cvm)] ])
    
    return(structure(list(
        call = this.call,
        cvm = cvm,
        cvsd = cvsd,
        lambda1 = lambda1,
        lambda2 = tn_init$lambda2,
        lambda1.min = lambda1.min,
        lambda1.1se = lambda1.1se
        ), class = "cv.tobitnet"))
}

plot.cv.tobitnet = function(x, ...){
    plot(x = log(x$lambda1), y = x$cvm, pch = 21, col = "red", bg = "red",
         xlab = expression( log(lambda [1]) ), ylab = "Mean-Squared Error", ...)
    arrows(log(x$lambda1), x$cvm-x$cvsd, log(x$lambda1), x$cvm+x$cvsd, 
           length=0.05, angle=90, code=3, col = "gray")
    abline(v = log(x$lambda1.min), lty = "dotted")
    abline(v = log(x$lambda1.1se), lty = "dotted")
}

plot.tobitnet = function(x, label = FALSE, ...){
    stopifnot( is.logical(label), length(label) == 1 )
    
    matplot(x = log(x$lambda1), y = t(x$beta), 
            xlab = expression( log(lambda[1]) ), ylab = "Coefficients",
            type = "l", lty = 1, ...)
    if(label == T){
        lmax = log(x$lambda1[length(x$lambda1)])
        coefmax = x$beta[,length(x$lambda1)]
        text(x = lmax, y = coefmax, labels = 1:length(coefmax), adj = c(2,0.5), cex = 0.5 )
    }
}
