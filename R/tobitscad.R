dscad = function(t, lambda, a){
  lambda*(t <= lambda) + (a*lambda - t)*(a*lambda > t)/(a - 1)*(t > lambda)
}

#3 step LLA with SCAD penalty
tobitscad = function(x, y, c = 0, a = 3, iter = 3, nlambda = 100, lambda.factor = ifelse(n/2 < nvars, 0.1, 0.05),
                     lambda = NULL, eps = 1e-7, standardize = TRUE, maxit = 1e6, early.stop = TRUE){
    this.call = match.call()
    
    n = nrow(x)
    nvars = ncol(x)
    varnames = colnames(x)
    
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
    }
    p = ncol(x)
    
    #c is numeric and scalar
    if(!is.numeric(c) | length(c) != 1 | !is.finite(c)) stop("c must be a finite scalar") 
    
    d = (y > c)
    
    #make sure a > 2 and iter >=2 (other argument checks are handled by tobitnet)
    stopifnot(is.numeric(a), is.finite(a), a > 2, length(a) == 1)
    stopifnot(is.numeric(iter), iter == round(iter), is.finite(iter), iter >= 2, length(iter) == 1)
    
    #Initial lasso fit
    tn1 = tobitnet(x, y, c = c, nlambda = nlambda, lambda.factor = lambda.factor,
                   lambda1 = lambda, lambda2 = 0, pf1 = rep(1, p), pf2 = rep(1, p),
                   eps = eps, standardize = standardize, maxit = maxit, early.stop = early.stop)
    
    #Pull initial deltas, lambda path (both needed to compute penalty factors)
    lambda_path = tn1$lambda1
    nlambda = length(lambda_path)
    delta = tn1$beta%*%diag(1/tn1$sigma, ncol = length(tn1$sigma))
    
    null_dev = tn1$nulldev
    dev = rep(0, nlambda)
    
    #LLA loop (need to complete for each lambda)
    for(i in 1:(iter-1)){
        gamma_vec = rep(0, nlambda)
        b0_vec = rep(0, nlambda)
        beta_mat = matrix(0, nrow = p, ncol = nlambda)
        
        delta_0_init = 0 
        delta_init = delta[,1]
        gamma_init = 1
        
        for(l in 1:nlambda){
            l1 = lambda_path[l]
            pfs = dscad(abs(delta[,l]), lambda = l1, a = a)/l1 
            
            #pass to tobitnet_innerC so pf1 can vary with lambda
            tn = tobitnet_innerC(xin = x, yin = y, cin = c, lambda1 = l1, lambda2 = 0, pf1 = pfs, pf2 = rep(1, p),
                                 delta_0_init = delta_0_init, delta_init = delta_init, gamma_init = gamma_init,
                                 eps = eps, standardize = standardize, maxit = maxit)
            if(!tn$KKT) warning("KKT conditions not satisfied. Try decreasing eps or increasing maxit.")
            
            if(i == iter - 1){
                #Compute deviance
                delta_pred = matrix(tn$beta*tn$gamma, nrow = p, ncol = 1)
                r_temp = x%*%delta_pred + rep((tn$b0 - c)*tn$gamma, n)
                dev[l] = -2*sum( logL1(y = y-c, d = d, r = r_temp, gamma = tn$gamma) )
            }
            
            gamma_vec[l] = tn$gamma
            b0_vec[l] = tn$b0
            beta_mat[,l] = tn$beta
            
            #Set previous solution as new inits (Note: tn$d0 and tn$delta are in the scale of the standardized predictors if standardize == T)
            delta_0_init = tn$d0
            delta_init = tn$delta
            gamma_init = tn$gamma
      }
      
      #overwrite delta for next LLA iteration
      delta = beta_mat%*%diag(gamma_vec, nrow = length(gamma_vec))
    }
    
    #Put in rows of 0s for betas corresponding to the constant columns
    beta_final = matrix(0, nrow = nvars, ncol = nlambda)
    beta_final[!constant_cols, ] = beta_mat
    rownames(beta_final) = varnames
    
    return(structure( list(
        call = this.call,
        sigma = 1/gamma_vec,
        b0 = b0_vec,
        beta = beta_final,
        c = c,
        lambda = lambda_path,
        dev = dev,
        nulldev = null_dev), class = "tobitscad")
        )
}

predict.tobitscad = function(object, newx, lambda = NULL, type = c("censored", "uncensored"), ...){
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
    
    if(!is.null(lambda)){
        object = update(object, lambda = lambda, c = c, ...)
    }
    
    beta_0 = matrix(object$b0, nrow = 1, ncol = length(object$lambda))
    
    r = newx%*%object$beta + matrix(1, nrow = n, ncol = 1)%*%beta_0
    if(type == 'censored'){
        r = pmax(r, c)
    }
    
    return(r)
}

cv.tobitscad = function(x, y, c = 0, a = 3, iter = 3, nlambda = 100, lambda.factor = ifelse(n/2 < nvars, 0.1, 0.05), 
                        lambda = NULL, nfolds = 10, early.stop = TRUE, type.measure = c("mse", "deviance", "mae"), ...){
    type.measure = match.arg(type.measure)
    this.call = match.call()
    n = nrow(x)
    nvars = ncol(x)
    p = nvars
    
    if( any( c(is.null(p), p <= 1) ) ) stop("x must be a matrix with 2 or more columns")
    
    #nfolds must be an integer between 2 and n (other argument checks covered by tobitnet and tobitscad)
    if(!is.numeric(nfolds) | nfolds != round(nfolds) | !(nfolds >= 2) | !(nfolds <= n) | !( is.finite(nfolds) ) | length(nfolds) != 1 ) stop("nfolds must be a single positive integer between 2 and the number of observations")
    
    tn_init = tobitnet(x = x, y = y, c = c, nlambda = nlambda, lambda.factor = lambda.factor, 
                       lambda1 = lambda, lambda2 = 0, pf1 = rep(1,p), pf2 = rep(1,p), early.stop = early.stop, ...)
    lambda = tn_init$lambda1
    nlambda = length(lambda)
    
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
        fold_size = length(foldlist[[i]])
        
        tscad = tobitscad(x = x[-foldlist[[i]], ], y = y[ -foldlist[[i]] ], c = c, a = a, iter = iter, lambda = lambda, early.stop = F, ...)
        
        if(type.measure == "mse"){
            preds = predict(tscad, newx = x[foldlist[[i]],], type = "censored") #n x nlambda
            r2 = (matrix(y[ foldlist[[i]] ], nrow = length(foldlist[[i]]), ncol = nlambda) - preds)^2
            err_mat[i,] = colMeans(r2)
        } else if(type.measure == "mae"){
            preds = predict(tscad, newx = x[foldlist[[i]],], type = "censored") #n x nlambda
            r_abs = abs(matrix(y[ foldlist[[i]] ], nrow = length(foldlist[[i]]), ncol = nlambda) - preds)
            err_mat[i,] = colMeans(r_abs)
        } else if(type.measure == "deviance"){
            delta_pred = tscad$beta/rep(tscad$sigma, each = nrow(tscad$beta))
            r_temp = as.matrix(x[ foldlist[[i]], ])%*%delta_pred + matrix( rep((tscad$b0 - c)/tscad$sigma, each = fold_size), nrow = fold_size, ncol = length(tscad$sigma))
            err_mat[i,] = vapply(1:ncol(r_temp), function(j) -2*sum(logL1(y = y[ foldlist[[i]] ]-c, d = (y[ foldlist[[i]] ] > c), r = r_temp[,j], gamma = 1/tscad$sigma[j])), numeric(1) ) 
        }
    }
    
    cvm = colMeans(err_mat) 
    cvvar = colMeans(err_mat*err_mat) - cvm^2 
    cvsd = sqrt(cvvar/nfolds)
    lambda.min = lambda[which.min(cvm)]
    lambda.1se = max(lambda[ cvm < min(cvm) + cvsd[which.min(cvm)] ])
    
    return(structure(list(
        call = this.call,
        cvm = cvm,
        cvsd = cvsd,
        lambda = lambda,
        lambda.min = lambda.min,
        lambda.1se = lambda.1se,
        type.measure = type.measure
    ), class = "cv.tobitscad"))
}

plot.cv.tobitscad = function(x, ...){
    if(x$type.measure == "mse"){ ylabel = "Mean-Squared Error" }
    else if(x$type.measure == "mae"){ ylabel = "Mean Absolute Error" }
    else if(x$type.measure == "deviance"){ ylabel = "Deviance" }
    plot(x = log(x$lambda), y = x$cvm, pch = 21, col = "red", bg = "red",
         xlab = expression( log(lambda [1]) ), ylab = ylabel, ...)
    arrows(log(x$lambda), x$cvm-x$cvsd, log(x$lambda), x$cvm+x$cvsd, 
           length=0.05, angle=90, code=3, col = "gray")
    abline(v = log(x$lambda.min), lty = "dotted")
    abline(v = log(x$lambda.1se), lty = "dotted")
}

plot.tobitscad = function(x, label = FALSE, ...){
    stopifnot( is.logical(label), length(label) == 1 )
    
    matplot(x = log(x$lambda), y = t(x$beta), 
            xlab = expression( log(lambda[1]) ), ylab = "Coefficients",
            type = "l", lty = 1, ...)
    if(label == T){
        lmax = log(x$lambda[length(x$lambda)])
        coefmax = x$beta[,length(x$lambda)]
        text(x = lmax, y = coefmax, labels = 1:length(coefmax), adj = c(2,0.5), cex = 0.5 )
    }
}