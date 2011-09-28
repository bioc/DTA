# Author: Sebastian Dümcke <duemcke@genzentrum.lmu.de>

tls = function(formula,D=NULL,T=NULL,precision=.Machine$double.eps){
    
    
    cl = match.call()
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    
    intercept = attr(mt,"intercept") == 1
    
    if(is.null(D)){
        D = diag(nrow(x))
    }
    if(is.null(T)){
        T = diag(ncol(x)+1)
    }
    
    coefs = wtls.householder(x,y,D,T,epsilon=precision)
    
    res = .residuals(x,y,coefs)
    r.sq = .r.squared(res,y,NULL)
    
    dn <- colnames(x)
    if (is.null(dn)) 
        dn <- paste("x", 1L:ncol(x), sep = "")
    
    colnames(coefs) = dn
    
    ret = list(coefficients=coefs[1,],residuals=res,call=cl,terms=mt,r.squared=r.sq)
    class(ret) = "lm"
    
    return(ret)
} 

##implement weigthed total least square regression according to Golub and Van Loan (1980) in SIAM J.Numer.Anal Vol 17 No.6
# A = design matrix/ data matrix
# b = observation
# D = diagonal weigth matrix for b
# T = diagonal weigth matrix for A (dim(T) = dim(A) + 1)
wtls.householder = function(A,b,D,T,epsilon=.Machine$double.eps){
    
    stopifnot(nrow(cbind(A,b)) > ncol(cbind(A,b)))
    ## check all(di > 0) && all(ti) > 0
    stopifnot(all(diag(D) > 0) && all(diag(T) > 0))
    
    #1 compute SVD of D[A|b]T
    s = svd(D %*% cbind(A,b) %*% T, nu = nrow(cbind(A,b)),nv=ncol(cbind(A,b)))
    
    #find index p of singular value that differs from the smallest one by less than a given amount epsilon 
    p = which(s$d <= rev(s$d)[1] + epsilon)[1]
    p = p:ncol(cbind(A,b))
    
    V = s$v[,p]
    if(length(p) == 1){
        
        y = V[1:(length(V)-1)]
        alpha = V[length(V)]
        
    }else{
        #compute Householder reflection and solve for y and alpha
        unit = c(rep(0,ncol(V)-1),1)
        c_n1 = as.matrix(V[nrow(V),])
        norm = as.matrix(norm(c_n1,type="F")*unit - c_n1) #/norm(norm(c_n1,type="F")*unit - c_n1,type="F")
        quot = as.vector(2 / (t(norm) %*% norm)) * (norm%*%t(norm))
        
        ##transpose of householder reflection
        Q1 = diag(1,nrow=nrow(quot),ncol=ncol(quot)) - quot
        
        mat = V %*% t(Q1)

#        stopifnot(all(mat[nrow(mat),max(1,ncol(mat)-1)] <= epsilon))
        if(!all(mat[nrow(mat),max(1,ncol(mat)-1)] <= epsilon)){
            warning(paste("Householder did not work, left lower values != 0:",mat[nrow(mat),max(1,ncol(mat)-1)]))
        }
        
        y = mat[,ncol(mat)][1:(nrow(mat)-1)]
        alpha = mat[,ncol(mat)][nrow(mat)]
    }
    
    if(alpha == 0){warning("TLS has no solution")}
    stopifnot(alpha != 0) #no solutions
    
    sol = as.matrix(-1/(alpha*rev(diag(T)[1]))) %*% diag(T)[1:(length(diag(T))-1)] * y 
    
    return(sol)
    
    
}

.residuals = function(x,y,coef){
    
    if(is.na(coef["(Intercept)"])){
        intcept=0
    }else{
        intcept = coef["(Intercept)"]
        coef = coef[-1]
    }
    
    #Hesse normal form of hyperplane
    factors = c(-1,coef)
    n = c(factors) / norm(as.matrix(factors),type="f")
    p = intcept / norm(as.matrix(factors),type="f")
    
    r = cbind(y,x) %*% n + p
    
    return(r)
}

.r.squared = function(x,y,coef){
    tss = sum((y - mean(y))^2)
#    yhat = apply(x%*%t(coef),1,sum)
#    ess = sum((yhat - mean(y))^2)
    rss = sum(x^2)
    
#    return(ess/tss)
    return(1 - (rss/tss))
}
# A = design matrix/ data matrix
# b = observation
# D = diagonal weigth matrix for b
# T = diabonal wigth matrix for A (dim(T) = dim(A) + 1)
wtls.solve = function(A,b,D,T){
    
    #check arguments
    if(!is.matrix(A)){A = as.matrix(A)}
    
    #1 compute SVD of DAT_1 (T_1 is the diagonal matrix with values t_1 ... t_n
    T1 = T[-nrow(T),-ncol(T)]
    s = svd(D %*% A %*% T1, nu = nrow(A),nv=ncol(A))
    
    Sigma_hat = diag(c(s$d,rep(0,ncol(A) - length(s$d))),ncol=ncol(A),nrow=nrow(A)) 
    sigma_hat = rev(s$d)[1]
    
    #2 solve LS: min ||D(AX-b)||_2^2
    ls = lm(b ~ I(D %*% A) + 0 ,x=TRUE)
#    rho = ls$x
#    rho = min(ls$residuals)
    xls = ls$coefficients
    dim(xls) = c(length(xls),1)
    rho = norm(D %*% (A %*% xls - b),type="F")^2 #is this right? is it really ||x||_2?
    
    #find sigma
    lambda = T[nrow(T),ncol(T)]
    C = t(s$u) %*% (D %*% b) #spaltenvector ?!
    
    psi = function(sigma){
        sigma^2 * ((1/lambda^2) + sum(C^2/(sigma_hat^2-sigma^2))) - rho
    }
    solv = uniroot(psi,c(0,sigma_hat))
    sigma = solv$root
    
    if(sigma  >= sigma_hat){
        warning("No solution to TLS")
        stop()
    }else{
        sol = T1 %*% s$v %*% solve(t(Sigma_hat) %*% Sigma_hat - sigma^2 * diag(ncol(Sigma_hat))) %*% t(Sigma_hat) %*% s$u %*% (D %*% b)
    }
    
    #calculate TLS sensitivity (higher = better, more stable)
    sigma_n1 = rev(svd(D %*% cbind(A,b) %*% T)$d)[1]
    sens = sigma_hat - sigma_n1
    
    print(paste("Sensitivity:",sens))
    return(list(coefficients=sol,sens=sens))
    
}

#tls = function(A,b){
#    
#    #check arguments
#    if(!is.matrix(A)){A = as.matrix(A)}
#    if(!is.matrix(b)){b = as.matrix(b)}
#    
#    s = svd(cbind(A,b), nu = nrow(cbind(A,b)),nv=ncol(cbind(A,b)))
#    
#    #take transose of s$v ???
#    V12 = t(s$v)[1:ncol(A),1:ncol(A)]
#    V22 = t(s$v)[(ncol(s$v) - ncol(b) + 1):ncol(s$v),(ncol(s$v) - ncol(b) + 1):ncol(s$v)]
#    
#    stopifnot(dim(V22) != c(ncol(b),ncol(b)))
#    
#    inv_V22 = solve(V22) #if error (V22 singular) problem has no solution
#
#    if(all(dim(inv_V22) == c(1,1))){
#        xtls = -V12 * inv_V22
#    }else{
#        xtls = -V12 %*% inv_V22
#    }
#    
#    return(xtls)
#}

tls.iter = function(A,y,delta=1e-10){
    
    if(!is.matrix(A)){
        A = as.matrix(A)
    }
    
    Nc = t(A) %*% cbind(A,y)
    N = Nc[,1:(ncol(Nc)-1)]
    c = Nc[,ncol(Nc)]
    
    k.start = 0
    eta.start = solve(N) %*% c
    
#    #1st iteration
#    k.i = t(y - A %*% eta.start) %*% (y - A %*% eta.start) / (1 + t(eta.start) %*% eta.start)
#    eta.i = solve(N - k.i %*% diag(ncol(A))) %*% c
    
    eta.prev = c(0,eta.start)
    eta.i = eta.start
    while(norm(as.matrix(eta.prev[2] - eta.prev[1]),type="F") >= delta){
        k.i = t(y - A %*% eta.i) %*% (y - A %*% eta.i) / (1 + t(eta.i) %*% eta.i)
        
        eta.i = solve(N - k.i %*% diag(ncol(A))) %*% c
        
        eta.prev[1] = eta.prev[2]
        eta.prev[2] = eta.i
    }
    
    var = k.i / (nrow(A)-ncol(A))
    
    return(list(coefficients=eta.i,var=var))
}
