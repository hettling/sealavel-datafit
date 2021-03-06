
## sum of squared residuals
ssr <- function(sim, data) {
    sum((data - sim)^2)
}

## log likelihood, given data, presiction, degrees of freedom.
## gives same result as stats::logLik, but can be used on 'raw' data
my.loglik <- function(dat, pred, df) {
    ssr <- sum((dat-pred)^2)
    N <- length(dat)
    0.5*(- N * (log(2 * pi) + 1 - log(N) +  log(ssr)))
}

## own implementation of AIC, could also use ?AIC
my.aic <- function(loglik, df){
    k <- 2
    -2 * loglik + k * df
}


## The following objective functions below are used to optimize using the ?optim function
## objective function for f1-f3
cost.3func <- function(p) {
    phase1 <- p['phase1']
    phase2 <- p['phase2']
    phase3 <- p['phase3']

    amp1 <- p['amp1']
    amp2 <- p['amp2']
    amp3 <- p['amp3']

    sim <- ff(amp1, t=T, f=f1, phase1)+ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3)
    sq <- ssr(sim, level)
    return (sq)
}

## objective function for f1-f4
cost.4func <- function(p) {
    phase1 <- p['phase1']
    phase2 <- p['phase2']
    phase3 <- p['phase3']
    phase4 <- p['phase4']

    amp1 <- p['amp1']
    amp2 <- p['amp2']
    amp3 <- p['amp3']
    amp4 <- p['amp4']

    sim <- ff(amp1, t=T, f=f1, phase1)+ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3) + ff(amp4, t=T, f=f4, phase4)
    sq <- ssr(sim, level)
    return (sq)
}

## objective function for f1-f5
cost.5func <- function(p) {
    phase1 <- p['phase1']
    phase2 <- p['phase2']
    phase3 <- p['phase3']
    phase4 <- p['phase4']
    phase5 <- p['phase5']

    amp1 <- p['amp1']
    amp2 <- p['amp2']
    amp3 <- p['amp3']
    amp4 <- p['amp4']
    amp5 <- p['amp5']

    sim <- ff(amp1, t=T, f=f1, phase1)+ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3) + ff(amp4, t=T, f=f4, phase4) + ff(amp5, t=T, f=f5, phase5)
    sq <- ssr(sim, level)

    return (sq)
}


## objective function for f1-f5 with frequency optimized
cost.5func.freq <- function(p) {

    phase1 <- p[1];#['phase1']
    phase2 <- p[2];#['phase2']
    phase3 <- p[3];#['phase3']
    phase4 <- p[4];#['phase4']
    phase5 <- p[5]#['phase5']

    if (max (c(phase1, phase2, phase3, phase4, phase5)) > 2*pi) {
        return (Inf)
    }

    amp1 <- p[6]#['amp1']
    amp2 <- p[7]#['amp2']
    amp3 <- p[8]#['amp3']
    amp4 <- p[9]#['amp4']
    amp5 <- p[10]#['amp5']

    freq5 <- p[11]#['freq5']

    cat("Freq5 : ", freq5, '\n');
    if (freq5 <= 0) {
        return (1e7)
    }

    sim <- ff(amp1, t=T, f=f1, phase1) + ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3) + ff(amp4, t=T, f=f4, phase4) + ff(amp5, t=T, f=freq5, phase5)
    sq <- ssr(sim, level)

    return (sq)
}

## objective function for f1-f5 with frequency optimized
sim.5func.freq <- function(p) {

    phase1 <- p[1];#['phase1']
    phase2 <- p[2];#['phase2']
    phase3 <- p[3];#['phase3']
    phase4 <- p[4];#['phase4']
    phase5 <- p[5]#['phase5']

    amp1 <- p[6]#['amp1']
    amp2 <- p[7]#['amp2']
    amp3 <- p[8]#['amp3']
    amp4 <- p[9]#['amp4']
    amp5 <- p[10]#['amp5']

    freq5 <- p[11]#['freq5']

    sim <- ff(amp1, t=T, f=f1, phase1) + ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3) + ff(amp4, t=T, f=f4, phase4) + ff(amp5, t=T, f=freq5, phase5)
    return (sim)
}


## grid search to find suitable starting points for the optimization with nls
## this is because nls is very sensitive to starting values
grid.nls <- function(formula, start, ...) {

    phases <- c('phase1', 'phase2', 'phase3', 'phase4', 'phase5')
    range.phases <- seq(0, 2*pi, length.out=3)
    amplitudes <- c('amp1', 'amp2', 'amp3', 'amp4', 'amp5')
    range.amplitudes <- seq(1, 100, length.out=3)
    range.f5 <- seq(1, 2000, length.out=3)

    ll <- list()
    for (p in phases) {
        ll[[p]] <- range.phases
    }
    for (a in amplitudes) {
        ll[[a]] <- range.amplitudes
    }
    ll[['freq5']] <- rev(range.f5)

    grid <- expand.grid(ll)
    grid <- grid[, names(start)]

    ## take random subset of grid, for performance reasons
    grid.size <- 20000
    grid <- grid[sample(1:nrow(grid))[1:grid.size],]

    ## set upper bound for phase
    upper <- start
    upper[names(upper)] <- 1e7
    upper[phases] <- 2*pi
    lower <- start
    lower[names(lower)] <- 0
    lower['freq5'] <- 1e-7

    results <- list()
    for (i in 1:nrow(grid)) {
        cat('starting parameter set # ', i, '\n')
        start[names(start)] <- grid[i,names(start)]
        start <- unlist(start)
        model <- try ( nls ( formula, start=start, lower=lower, upper=upper, algorithm='port', control=nls.control(maxiter = 500,printEval=FALSE)))
        if (!inherits(model, "try-error")) {
            cat('optimization succesful\n')
            ##return(model);
            results[[length(results)+1]] <- model
        }
    }
    ## look at best performing model
    ssrs <- sapply(results, function(f)f$m$deviance())
    best.idx <- which(ssrs==min(ssrs))[1]
    best <- results[[best.idx]]
    cat ("SSR of best model (own grid search)  : ", best$m$deviance(), "\n")
    ## do another fit:
    reest <- try (nls(formula, start=best$m$getPars(), algorithm='port', lower=lower, upper=upper))
    if (!inherits(reest, "try-error")) {
        cat ("SSR of re-fitted best model (own grid search)  : ", reest$m$deviance(), "\n")
        return (reest)
    } else {
        return (best)
    }
    #require('nls2')
    ## make fit using grid nls2
    #model.2 <- nls2 ( formula, start=grid, lower=lower, upper=upper, algorithm='grid-search', control=nls.control(maxiter = 500,printEval=FALSE))
    #m.2 <- nls2(formula, start=model.2, algorithm='port', lower=lower, upper=upper)
    #cat ("SSR of best model (nls2)  : ", m.2$m$deviance(), "\n")

}

## fit models using nls, call the grid search to get good starting values
fit.nls <- function(data) {
    level <<- data
    cat("fitting using nls\n")

    ## fit with four orbital cycles
    fit.4func <- grid.nls(level~ff(amp1, t=T, f=f1, phase1)+ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3)+ff(amp4, t=T, f=f4, phase4),
                          start=list(phase1=0, phase2=0, phase3=0, phase4=0, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4))
    cat("fitted 4 func\n")

    ## fit with four orbital cycles, third order, fixed frequency
    fit.5func <- grid.nls(level~ff(amp1, t=T, f=f1, phase1)+ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3)+ff(amp4, t=T, f=f4, phase4)+ff(amp5, t=T, f=f5, phase5),
                          start=list(phase1=0, phase2=0, phase3=0, phase4=0, phase5=0, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4, amp5=amp5))
    cat("fitted 5 func\n")

    ## fit with four orbital cycles, third order including frequency
    fit.5func.freq <- grid.nls(level~ff(amp1, t=T, f=f1, phase1)+ff(amp2, t=T, f=f2, phase2)+ff(amp3, t=T, f=f3, phase3)+ff(amp4, t=T, f=f4, phase4)+ff(amp5, t=T, f=freq5, phase5),
                               start=list(phase1=0, phase2=0, phase3=0, phase4=0, phase5=0, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4, amp5=amp5, freq5=f5))
    cat("fitted 5 func with frequency \n")

    res <- list(fit.4func, fit.5func, fit.5func.freq)
    names (res) <- c("funcs4", "funcs5", "func5freq")
    return (res)
}


#    out <- adaptMetropGibbs(cost.mcmc, log(start.f5.freq), batch=10, batch.length=50000, report=1, accept.rate=(rep(0.25, length(start.f5.freq))), verbose=T)
## omit burnin
#    n <- out$batch*out$batch.length
#    burnin <- n / 10
#    ens <- exp(out$p.theta.samples)[burnin:n,]

## do a mcmc for model with tectonic cyclicity
mcmc.f5.freq <- function() {
    start.f5.freq <- c(phase1=1, phase2=2, phase3=1, phase4=2, phase5=1, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4, amp5=amp5, freq5=f5)
    cost.mcmc <- function(p) {
        cost <- cost.5func.freq(exp(p))
        return ((-cost))
    }
    ## we are in logspace, therefore have to scale to small steps to get good acceptance
    length <- 20000
    ## m <- metrop(cost.mcmc, log(start.f5.freq), length, scale=0.007)
    out <- adaptMetropGibbs(cost.mcmc, log(start.f5.freq), batch=10, batch.length=5000, report=1, accept.rate=(rep(0.25, length(start.f5.freq))), verbose=T)

    post <- exp(out$p.theta.samples)

    ssr <- apply(post, 1, cost.5func.freq)
    min.ssr <- min(ssr)
    cat("Min SSR for MCMC : ", min.ssr, "\n")

    ## discard first 10% as burn-in
    burnin <- (dim(post)[1] / 10)+ 1
    post <- post[burnin:dim(post)[1],]

    ## calculate autocorrelation times for thinning
    acf.mat <- apply(post, 2, function(x)acf(x, lag.max=10000, plot=FALSE, n.used=nrow(post))$acf)
    corrtimes <- apply( acf.mat , 2, function(x){ which(x < 1/exp(1))[1] })
    max.cor <- max(corrtimes)

    ## thinning
    ppost <- post[seq(1, dim(post)[1], max.cor),]

    return (post)
}

## fit models using the optim function
fit.optim <- function() {

    start.f4 <- c(phase1=0, phase2=0, phase3=0, phase4=0, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4)
    opt.f4 <- optim(start.f4, cost.4func, control=list(maxit=1e7), lower=0, method='L-BFGS-B')

    start.f5 <- c(phase1=0, phase2=0, phase3=0, phase4=0, phase5=0, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4, amp5=amp5)
    opt.f5 <- optim(start.f5, cost.5func, control=list(maxit=1e7), lower=0, method='L-BFGS-B')

    start.f5.freq <- c(phase1=1, phase2=2, phase3=1, phase4=2, phase5=1, amp1=amp1, amp2=amp2, amp3=amp3, amp4=amp4, amp5=amp5, freq5=f5)
    lower <- c(rep(0, length(start.f5.freq)-1), 1)
    upper <- c(rep(2*pi, 5), rep(10000, 6) )
    opt.f5.freq <- optim(start.f5.freq, cost.5func.freq, control=list(maxit=1e7), lower=lower, upper=upper, method='L-BFGS-B')

    res <- list(opt.f4, opt.f5, opt.f5.freq)
    names (res) <- c("funcs4", "funcs5", "funcs5.freq")
    return(res)
}


## predict sea level using optimized model. These functions are to predict when optimizing with optim
predict.3func <- function(par) {
    return( ff(par['amp1'], t=T, f=f1, par['phase1'])
           + ff(par['amp2'], t=T, f=f2, par['phase2'])
           +  ff(par['amp3'], t=T, f=f3, par['phase3']))
}

predict.4func <- function(par) {
    return(ff(par['amp1'], t=T, f=f1, par['phase1'])
           + ff(par['amp2'], t=T, f=f2, par['phase2'])
           +  ff(par['amp3'], t=T, f=f3, par['phase3'])
           +  ff(par['amp4'], t=T, f=f4, par['phase4']))
}

predict.5func <- function(par) {
    return(ff(par['amp1'], t=T, f=f1, par['phase1'])
           + ff(par['amp2'], t=T, f=f2, par['phase2'])
           +  ff(par['amp3'], t=T, f=f3, par['phase3'])
           +  ff(par['amp4'], t=T, f=f4, par['phase4'])
           +  ff(par['amp5'], t=T, f=f5, par['phase5']))
}

predict.5func.freq <- function(par) {
    return(ff(par['amp1'], t=T, f=f1, par['phase1'])
           + ff(par['amp2'], t=T, f=f2, par['phase2'])
           +  ff(par['amp3'], t=T, f=f3, par['phase3'])
           +  ff(par['amp4'], t=T, f=f4, par['phase4'])
           +  ff(par['amp5'], t=T, f=par['freq5'], par['phase5']))
}

