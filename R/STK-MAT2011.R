require(ggplot2)
# Formatting functions ####
date_time <- function(dates, times, start_date) {
  first <- as.POSIXct("1970-01-01 00:00:00")
  start_date <- as.POSIXct("2016-01-01 00:00:00")
  orig <- difftime(start_date,first, units = "secs")
  dts <- paste(dates, times)
  dts1 <- as.POSIXct(dts, format = "%d.%m.%Y %H:%M:%S")
  dts2 <- dts1 - orig
  dts3 <- as.numeric (dts2)
  dts4 <- dts3/(60*60)
}

# Plotting functions ####
cumfreqs <- function(time) {
  bin.width <- 1
  min <- min(time)
  max <- max(time)
  breaks <- seq(min, max+1, by = bin.width) # specify end points of bins
  cuts <- cut(time, breaks, right=FALSE) # assign each interval with a bin
  freqs <- table(cuts)
  cumFreq <- cumsum(freqs)
  DG_x <- c(min(time) : max(time))
  DG_y <- as.vector(cumFreq)
  tag <- c('time')
  DG_legend <- as.factor(rep(tag, times=length(DG_y)))
  DG_data <- data.frame(DG_x,DG_y,DG_legend)
  ggplot(data = DG_data, aes(x=DG_x, y=DG_y)) + geom_line(aes(color=DG_legend)) + xlab("Time") + ylab("Cumulative Number of Failures")
}


# Plots inter-failure times against observation times
ift_cf <- function(name, time, b0, b1, gam, delt) {
  plot(time, 1:length(time), pch='', ylab = "Index", xlab = "Time in hrs", main=name)
  i <- seq(1, length(time), 1)
  segments(time[i], i, time[i+1], i)
  lines(time, mu_1(time, b0, b1), col='blue')
  lines(time, mu_2(time, gam, delt), col='magenta')
}
# Estimation functions ####

# Fitting the model, HPP
mle.hpp <- function(time) {
  lambda_hat <- mean(time)
  return(lambda_hat)
}

loglik.hpp <- function(time, lambda) {as.numeric(exp(-lambda*length(time))*lambda^sum(time))}

# Fitting the model, NHPP

#Log-linear
mle.ll <- function(time) {
  n <- length(time)
  time.end <- time[length(time)]
  ll.lik <- function(b_1)
  {sum(time) + n/b_1 - (n*time.end)/(1-exp(-b_1*time.end))}
  b_1 <- uniroot(ll.lik, interval=c(-1,1))$root
  b_0 <- log(b_1*n/(exp(b_1*time.end)-1))
  return(c(b_0, b_1))
}

lambda_ll <- function(time, b_0, b_1) {exp(b_0+b_1*time)}
mu_ll <- function(time, b_0, b_1) {exp(b_0)/b_1*(exp(b_1*time)-1)}

# Power law
mle.pl <- function(time) {
  n <- length(time)
  time.end <- time[length(time)]
  delt <- n/(n*log(time.end)-sum(log(time)))
  gam <- n/time.end^delt
  return(c(gam, delt))
}

lambda_pl <- function(time, gam, delt) {gam*delt*time^(delt-1)}
mu_pl <- function(time, gam, delt) {gam*time^delt}

#Max loglik
loglik.ll <- function(time, b_0, b_1) {length(time)*b_0 + b_1*sum(time)-length(time)}
loglik.pl <- function(time, gam, delt) {(delt-1)*sum(log(time)) + length(time)*log(gam)+length(time)*log(delt) - length(time)}



## The Laplace Test
laplace_test <- function(times) {
  n <- length(times)
  t.e <- times[length(times)]
  z <- (sqrt(12*n)*sum(times-t.e/2))/(n*t.e)
  return(z)
}

# Military handbook test
military_test <- function(times) {
  n <- length(times)
  t.e <- times[length(times)]
  z <- 2*sum(log(t.e/times))
  chi <- pchisq(z, 2*n)
  return(c(z, chi))
}


