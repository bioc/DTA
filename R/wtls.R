

### implement weigthed total least square regression according to Golub and Van Loan (1980) in SIAM J.Numer.Anal Vol 17 No.6 ###


tls = function(formula,
		D=NULL,
		T=NULL,
		precision=.Machine$double.eps
)
{
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
    if(is.null(D)){D = diag(nrow(x))}
    if(is.null(T)){T = diag(ncol(x)+1)}
    coefs = wtls.householder(x,y,D,T,epsilon=precision)
    res = .residuals(x,y,coefs)
    r.sq = .r.squared(res,y,NULL)
    dn <- colnames(x)
    if (is.null(dn)) dn <- paste("x", 1L:ncol(x), sep = "")
    colnames(coefs) = dn
    ret = list(coefficients=coefs[1,],residuals=res,call=cl,terms=mt,r.squared=r.sq)
    class(ret) = "lm"
    return(ret)
} 


### wtls.householder function ###


wtls.householder = function(A,
		b,
		D,
		T,
		epsilon=.Machine$double.eps
)
{
    stopifnot(nrow(cbind(A,b)) > ncol(cbind(A,b)))
    stopifnot(all(diag(D) > 0) && all(diag(T) > 0))
    s = svd(D %*% cbind(A,b) %*% T, nu = nrow(cbind(A,b)),nv=ncol(cbind(A,b)))
    p = which(s$d <= rev(s$d)[1] + epsilon)[1]
    p = p:ncol(cbind(A,b))
    V = s$v[,p]
    if(length(p) == 1){
        y = V[1:(length(V)-1)]
        alpha = V[length(V)]
    }else{
        unit = c(rep(0,ncol(V)-1),1)
        c_n1 = as.matrix(V[nrow(V),])
        norm = as.matrix(norm(c_n1,type="F")*unit - c_n1)
        quot = as.vector(2 / (t(norm) %*% norm)) * (norm%*%t(norm))
        Q1 = diag(1,nrow=nrow(quot),ncol=ncol(quot)) - quot
        mat = V %*% t(Q1)
	if(!all(mat[nrow(mat),max(1,ncol(mat)-1)] <= epsilon)){
            warning(paste("Householder did not work, left lower values != 0:",mat[nrow(mat),max(1,ncol(mat)-1)]))
        }
        y = mat[,ncol(mat)][1:(nrow(mat)-1)]
        alpha = mat[,ncol(mat)][nrow(mat)]
    }
    if(alpha == 0){warning("TLS has no solution")}
    stopifnot(alpha != 0)
    sol = as.matrix(-1/(alpha*rev(diag(T)[1]))) %*% diag(T)[1:(length(diag(T))-1)] * y 
    return(sol)
}


### .residuals function ###


.residuals = function(x,y,coef)
{
    if(is.na(coef["(Intercept)"])){
        intcept=0
    }else{
        intcept = coef["(Intercept)"]
        coef = coef[-1]
    }
    factors = c(-1,coef)
    n = c(factors) / norm(as.matrix(factors),type="f")
    p = intcept / norm(as.matrix(factors),type="f")
    r = cbind(y,x) %*% n + p
    return(r)
}


### .r.squared function ###


.r.squared = function(x,y,coef)
{
    tss = sum((y - mean(y))^2)
    rss = sum(x^2)
    return(1 - (rss/tss))
}


### wtls.solve function ###


wtls.solve = function(A,b,D,T)
{
   	if(!is.matrix(A)){A = as.matrix(A)}
    T1 = T[-nrow(T),-ncol(T)]
    s = svd(D %*% A %*% T1, nu = nrow(A),nv=ncol(A))
    Sigma_hat = diag(c(s$d,rep(0,ncol(A) - length(s$d))),ncol=ncol(A),nrow=nrow(A)) 
    sigma_hat = rev(s$d)[1]
    ls = lm(b ~ I(D %*% A) + 0 ,x=TRUE)
    xls = ls$coefficients
    dim(xls) = c(length(xls),1)
    rho = norm(D %*% (A %*% xls - b),type="F")^2
    lambda = T[nrow(T),ncol(T)]
    C = t(s$u) %*% (D %*% b)
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
    sigma_n1 = rev(svd(D %*% cbind(A,b) %*% T)$d)[1]
    sens = sigma_hat - sigma_n1
    print(paste("Sensitivity:",sens))
    return(list(coefficients=sol,sens=sens))   
}


### tls.iter function ###


tls.iter = function(A,y,delta=1e-10)
{
    if(!is.matrix(A)){
        A = as.matrix(A)
    }
    Nc = t(A) %*% cbind(A,y)
    N = Nc[,1:(ncol(Nc)-1)]
    c = Nc[,ncol(Nc)]
    k.start = 0
    eta.start = solve(N) %*% c
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



