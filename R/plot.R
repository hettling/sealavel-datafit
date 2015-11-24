pdf(plotfile, family="Times", width=3.4, height=3.4, pointsize=7)
par(mar=c(4.1, 4.1, 1.5, 1.1))
plot(level~T, ylim=c(-80, 120), xlab="time (k years)", ylab="eustatic sea level (m)", type='n', cex.lab=1.3)

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


# plot data
points(level~T, xlab="time (k years)", ylab="eustatic sea level (m)", pch=18)

# plot fit with four orbital cycles
lines(T, sim.nls$funcs4, col='blue', lwd=2)

# plot fit with four orbital and third order cycle
lines(T, sim.nls$func5freq, col='red', lwd=2)

legend('bottomleft', col=c('blue', 'red'), legend=c('orbital cycles', 'orbital cycles + third order cycle'), box.lwd=0, lwd=2, box.col='white')
dev.off()
