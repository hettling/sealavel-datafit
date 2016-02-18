

require('lmtest')
set.seed(111)

datafile <- 'sealevel-comp-max.txt'
#datafile <- 'sealevel-comp-min.txt'
#datafile <- 'sealevel-comp-max.txt'
plotfile <- paste0(substr(datafile, 0, nchar(datafile)-4), '.pdf')
restable <- paste0(substr(datafile, 0, nchar(datafile)-4), '.tsv')

dat <- read.table(paste0('../data/', datafile), header=T)

#dat <- dat[which(dat$interval %in% c("I", "II", "III")),]

source('funcs.R')

### cycle amplitutes that can be used as starting values for optimization
amp1 <- 1#10
amp2 <- 1#40
amp3 <- 1#50
amp4 <- 1#50
amp5 <- 1000

## cycle frequencies
f1 <- 23
f2 <- 41
f3 <-100
f4 <- 412
f5 <- 2000 #initial value of third order, gets estimated

## This is our sinusoidal function
## For multiple cycles, multiple instances of this function (with different parameters) will be summed
ff <- function(amp, t, f, phase) {
    amp * sin(2*pi * 1/f * t + phase)
}

## get data, units : kyears and meters, respectively
#interval <- dat[,'interval']
T <- dat[,'time']/1000
offset <- min(T) 
T <- T - offset 
level <- dat[,'level']


## Results below are named lists of length 5 consisting the following models:
##  1)  funcs4: model with four orbital cycles
##  2)  funcs5: model with four orbital cycles and third order (fixed frequency of third order)
##  3)  functs5.freq: model with four orbital cycles, third order and optimized third order frequency

## fits for all models with nls
res.nls <- fit.nls()

## fits for all models with optim
res.opt <- fit.optim()

## predict sea level with all models
sim.nls <- lapply(res.nls, predict)
sim.opt <- list(predict.4func(res.opt[[1]]$par), predict.5func(res.opt[[2]]$par), predict.5func.freq(res.opt[[3]]$par))

# get sum of squared residuals for all models
ssr.nls <- lapply(sim.nls, ssr, level)
ssr.opt <- lapply(sim.opt, ssr, level)

## get degrees of freedom for all models
degrees <- list(funcs4=8, funcs5=10, func5freq=11)# <- lapply(res.nls, function(x){summary(x)$df[1]})

## calculate log likelihood
loglik.nls <- mapply(my.loglik, list(level, level, level), sim.nls, degrees)
loglik.opt <- mapply(my.loglik, list(level, level, level), sim.opt, degrees)

## calculate AIC
aic.nls <- mapply(my.aic, loglik.nls, degrees)
aic.opt <- mapply(my.aic, loglik.opt, degrees)

## likelihood ratio tests between various models
t4.5 <- lrtest(res.nls[['funcs4']], res.nls[['funcs5']])
t4.5.fr <- lrtest(res.nls[['funcs4']], res.nls[['func5freq']])
t5.5.fr <- lrtest(res.nls[['funcs5']], res.nls[['func5freq']])

## change to negative times again
T <- T + offset

## make plot for figure 20
source('plot.R')

## print values needed for paper
pars <- res.nls$func5freq$m$getPars()

## frequency of fifth cycle
#cat("Tectonic cyclicity : ", pars['freq5'], " years\n")

## results of likelihood ratio test
#print(t4.5.fr)

## AIC
#cat("AIC simpler model: ", aic.nls[2], " tectonic model : ", aic.nls[3], "\n");

## get p value for likelihood rato test
pval <- t4.5.fr["Pr(>Chisq)"][[1]][2]
loglik.4func <- round(t4.5.fr$LogLik[1], 3)
loglik.5func.freq <- round(t4.5.fr$LogLik[2], 3)

aic.4func <- round(aic.nls[2], 3)
aic.5func.freq <- round(aic.nls[3], 3)

## amplitudes
##cat("Amplitudes : ", round(pars[paste0('amp', 1:5)],1), "\n")

## write results to table
results <- c(tectonic_cyclicity=unname(pars['freq5']), round(pars[paste0('amp', 1:5)],1), round(pars[paste0('phase', 1:5)],1), ssr=round(ssr.nls[[3]],2), aic_4func=aic.4func, aic_5func=aic.5func.freq, loglik_4func=loglik.4func, loglik_5func=loglik.5func.freq, pvalue=pval)
write.table(data.frame(t(results)), file=restable, quote=F, sep="\t", row.names=F)


