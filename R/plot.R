pdf(plotfile, family="Times", width=3.4, height=3.4, pointsize=7)
par(mar=c(4.1, 4.1, 1.5, 1.1))
plot(level_avg~T, ylim=c(-80, 120), xlab="time (k years)", ylab="eustatic sea level (m)", type='n', cex.lab=1.3)

color <- 'grey80'

# plot intervals as shaded regions
for (iv in unique(dat$interval)) {
    idx <- which(dat$interval==iv)

    # avoid spaces; interval ends where new interval starts!!
    if (iv != dat$interval[nrow(dat)]) {
        idx <- c(idx, tail(idx, 1)+1)
    }
    t <- dat$time[idx]/1000

    rect(min(t), -60, max(t), 110, col=color, lty=0)
    color <- ifelse(color=='grey80', 'grey90', 'grey80')
    text((min(t)+max(t))/2, 115, labels=iv, cex=0.7)
}


## plot data
points(level_avg~T, xlab="time (k years)", ylab="eustatic sea level (m)", pch=18)
## error bars
arrows(T, level_avg - level_sds, T, level_avg + level_sds, length=0.02, angle=90, code=3, lwd=0.2)

## plot fit with four orbital cycles
lines(T, sim.nls.avg$funcs4, col='blue', lwd=2)
## upper and lower bound
lines(T, sim.nls.max$funcs4, col='blue', lwd=1, lty=2)
lines(T, sim.nls.min$funcs4, col='blue', lwd=1, lty=2)

## plot fit with four orbital and third order cycle
lines(T, sim.nls.avg$func5freq, col='red', lwd=2)
## upper and lower bound
lines(T, sim.nls.max$func5freq, col='red', lwd=1, lty=2)
lines(T, sim.nls.min$func5freq, col='red', lwd=1, lty=2)

legend('bottomleft', col=c('blue', 'red'), legend=c('orbital cycles', 'orbital cycles + third order cycle'), box.lwd=0, lwd=2, box.col='white')
dev.off()


## plot individual cycles

##pars <- res$m$getPars()
#plot(T, dat$level)
#lines(T, sim.nls[[3]], col='black')

#res <- res.nls[[3]]
#sim <- sim.nls[[3]]

##time <- seq(min(T), max(T), 0.1)
#time <- seq(0, 800, 0.1)
#c1 <- ff(pars['amp1'], t=time, f=f1, phase=pars['phase1'])
#lines(time, c1, col='blue', type='l')

#c2 <- ff(pars['amp2'], t=time, f=f2, phase=pars['phase2'])
#lines(time, c2, col='red')

#c3 <- ff(pars['amp3'], t=time, f=f3, phase=pars['phase3'])
#lines(time, c3, col='green')

#c4 <- ff(pars['amp4'], t=time, f=f4, phase=pars['phase4'])
#lines(time, c4, col='yellow')

#c5 <- ff(pars['amp5'], t=time, f=pars['freq5'], phase=pars['phase5'])
#lines(time, c5, col='brown')



