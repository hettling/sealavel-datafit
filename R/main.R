require('lmtest')
require('mcmc')
set.seed(111)

datafile <- 'sealevel-comp-min.txt'
#datafile <- 'sealevel-comp-avg.txt'
#datafile <- 'sealevel-comp-max.txt'
plotfile <- paste0(substr(datafile, 0, nchar(datafile)-4), '.pdf')
restable <- paste0(substr(datafile, 0, nchar(datafile)-4), '.tsv')

dat <- read.table(paste0('../data/', datafile), header=T, stringsAsFactors=F)

#dat <- dat[which(dat$interval %in% c("I", "II", "III")),]

source('funcs.R')

### cycle amplitutes that can be used as starting values for optimization
amp1 <- 1#10
amp2 <- 1#40
amp3 <- 1#50
amp4 <- 1#50
amp5 <- 10

## cycle frequencies
# first two values from Berger & Loutre (1994)
f1 <- 19 # 23 # axial precission (trend of earth axis direction to fix stars)
f2 <- 38.5 #41 # obliquity (axial tilt)
f3 <-100 # eccentricity 1 (shape of earths orbit around sun)
f4 <- 412 # eccentricity 2
f5 <- 50#2000 #initial value of third order, gets estimated

## This is our sinusoidal function
## For multiple cycles, multiple instances of this function (with different parameters) will be summed
ff <- function(amp, t, f, phase) {
    amp * sin(2*pi * 1/f * t + phase)
}

## get data, units : kyears and meters, respectively
#interval <- dat[,'interval']
T <- as.vector(dat[,'time'])/1000
offset <- min(T)
T <- T - offset
level_min <- dat[,'level_min']
level_max <- dat[,'level_max']
level_avg <- dat[,'level_avg']

level_sds<- apply(cbind(level_min, level_max), 1, sd)


level <- level_avg

## Results below are named lists of length 5 consisting the following models:
##  1)  funcs4: model with four orbital cycles
##  2)  funcs5: model with four orbital cycles and third order (fixed frequency of third order)
##  3)  functs5.freq: model with four orbital cycles, third order and optimized third order frequency

## fits for all models with nls
res.nls.avg <- fit.nls(level_avg)
res.nls.min <- fit.nls(level_min)
res.nls.max <- fit.nls(level_max)

## fits for all models with optim
##res.opt <- fit.optim()

## predict sea level with all models
sim.nls.avg <- lapply(res.nls.avg, predict)
sim.nls.min <- lapply(res.nls.min, predict)
sim.nls.max <- lapply(res.nls.max, predict)

##sim.opt <- list(predict.4func(res.opt[[1]]$par), predict.5func(res.opt[[2]]$par), predict.5func.freq(res.opt[[3]]$par))

# get sum of squared residuals for all models
ssr.nls <- lapply(sim.nls.avg, ssr, level_avg)
##ssr.opt <- lapply(sim.opt, ssr, level)

## get degrees of freedom for all models
degrees <- list(funcs4=8, funcs5=10, func5freq=11)# <- lapply(res.nls, function(x){summary(x)$df[1]})

## calculate log likelihood
loglik.nls <- mapply(my.loglik, list(level, level, level), sim.nls.avg, degrees)
##loglik.opt <- mapply(my.loglik, list(level, level, level), sim.opt, degrees)

## calculate AIC
aic.nls <- mapply(my.aic, loglik.nls, degrees)
##aic.opt <- mapply(my.aic, loglik.opt, degrees)

## likelihood ratio tests between various models
t4.5 <- lrtest(res.nls.avg[['funcs4']], res.nls.avg[['funcs5']])
t4.5.fr <- lrtest(res.nls.avg[['funcs4']], res.nls.avg[['func5freq']])
t5.5.fr <- lrtest(res.nls.avg[['funcs5']], res.nls.avg[['func5freq']])

## change to negative times again
T <- T + offset

## make plot for figure 20
source('plot.R')

## print values needed for paper
pars.avg <- round(res.nls.avg$func5freq$m$getPars(), 2)
pars.min <- round(res.nls.min$func5freq$m$getPars(), 2)
pars.max <- round(res.nls.max$func5freq$m$getPars(), 2)

## get p value for likelihood rato test
pval <- t4.5.fr["Pr(>Chisq)"][[1]][2]
loglik.4func <- round(t4.5.fr$LogLik[1], 3)
loglik.5func.freq <- round(t4.5.fr$LogLik[2], 3)

aic.4func <- round(aic.nls[2], 3)
aic.5func.freq <- round(aic.nls[3], 3)


## Generate table

pars.total <- vector()

## Generate ranges
for (n in names(pars.avg)) {
    str <- paste(pars.avg[n], '(', min(pars.min[n], pars.avg[n], pars.max[n]), '-', max(pars.min[n], pars.avg[n], pars.max[n]), ')')
    pars.total[n] <- str
}

## write results to table
results <- c(
    tectonic_cyclicity=unname(pars.total['freq5']),
    pars.total[paste0('amp', 1:5)],
    pars.total[paste0('phase', 1:5)],
    ssr=round(ssr.nls[[3]],2),
    aic_4func=aic.4func,
    aic_5func=aic.5func.freq,
    loglik_4func=loglik.4func,
    loglik_5func=loglik.5func.freq,
    pvalue=pval)
write.table(data.frame(t(results)), file=restable, quote=F, sep="\t", row.names=F)


