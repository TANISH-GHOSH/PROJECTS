# ==============================================================================
# SECTION 1: Monte Carlo Principle Demonstration
# Description: Demonstrates the Law of Large Numbers by simulating standard 
# normal variables and verifying the sample mean and variance.
# ==============================================================================

set.seed(123)
N <- 10000
samples <- rnorm(N, mean = 0, sd = 1)

mean(samples)   # should be close to 0
var(samples)    # should be close to 1


# ==============================================================================
# SECTION 2: Brownian Motion Simulation
# Description: Defines a function to simulate general Brownian Motion with drift 
# and volatility, and plots four different parameter configurations.
# ==============================================================================

set.seed(123)

sim.Brownian.motion = function(X0, t, mu, sigma){
  n = length(t)
  X = rep(NA, n)
  Z = rnorm(n, mean=0, sd=1)
  X[1] = X0 + sigma*sqrt(t[1])*Z[1] + mu*t[1]
  for(i in 2:n){
    X[i] = X[i-1] + sigma*sqrt(t[i]-t[i-1])*Z[i] + mu*(t[i]-t[i-1])
  }
  return(X)
}

t = seq(0, 1, 0.001)

par(mfrow = c(2,2))

# Top Left: X0=0, mu=0, sigma=1
X1 = sim.Brownian.motion(X0=0, t=t, mu=0, sigma=1)
plot(ts(X1, start=0, frequency=1000), ylab="",
     main=expression(X[0]==0~","~mu==0~","~sigma==1))

# Top Right: X0=0, mu=2, sigma=1
X2 = sim.Brownian.motion(X0=0, t=t, mu=2, sigma=1)
plot(ts(X2, start=0, frequency=1000), ylab="",
     main=expression(X[0]==0~","~mu==2~","~sigma==1))

# Bottom Left: X0=0, mu=0, sigma=2
X3 = sim.Brownian.motion(X0=0, t=t, mu=0, sigma=2)
plot(ts(X3, start=0, frequency=1000), ylab="",
     main=expression(X[0]==0~","~mu==0~","~sigma==2))

# Bottom Right: X0=2, mu=0, sigma=1
X4 = sim.Brownian.motion(X0=2, t=t, mu=0, sigma=1)
plot(ts(X4, start=0, frequency=1000), ylab="",
     main=expression(X[0]==2~","~mu==0~","~sigma==1))

par(mfrow = c(1,1))


# ==============================================================================
# SECTION 3: Brownian Bridge Construction
# Description: Simulates a standard Brownian path, pins the terminal value, 
# and constructs a Brownian bridge interpolating between endpoints.
# ==============================================================================

set.seed(123)

T       <- 1
n_steps <- 100
dt      <- T / n_steps
t_grid  <- seq(0, T, length.out = n_steps)

# Simulate a standard Brownian path
Z <- rnorm(n_steps - 1)
W <- c(0, cumsum(sqrt(dt) * Z))

# Fix terminal value and construct the bridge
W_T <- W[n_steps]
BB  <- W - (t_grid / T) * W_T   # B(t) = W(t) - (t/T)*W(T)

plot(t_grid, BB, type = "l",
     main = "Brownian Bridge (pinned at zero)",
     xlab = "Time", ylab = "B(t)")
points(c(0, T), c(0, 0), pch = 19, col = "red")
legend("topright", legend = "Fixed endpoints", pch = 19, col = "red")


# ==============================================================================
# SECTION 4: Principal Components Construction
# Description: Uses PCA on the covariance matrix of Brownian motion increments 
# to simulate paths with full and truncated eigenvectors (smoothing effect).
# ==============================================================================

set.seed(896)
sig = matrix(0,1000,1000)
for(i in 1:1000){
  for(j in 1:1000){
    sig[i,j] = min(i,j)/1000
  }
}

u = eigen(sig)
v = cumsum(u$values)
v = v/v[1000]
z = rnorm(1000,0,1)
zz = z*sqrt(u$values)
ful = u$vectors%*%zz

app20  = u$vectors[,1:20]  %*% zz[1:20]
app100 = u$vectors[,1:100] %*% zz[1:100]

plot(ts(ful, start=0, frequency=1000), ylab="", lwd=2)
lines(ts(app20,  start=0, frequency=1000), col='red',  lwd=2)
lines(ts(app100, start=0, frequency=1000), col='blue', lwd=2)
legend("topright",
       legend = c("Full (1000 components)", 
                  "Approximation (100 components)",
                  "Approximation (20 components)"),
       col    = c("black", "blue", "red"),
       lwd    = 2)


# ==============================================================================
# SECTION 5: Geometric Brownian Motion (GBM) - Single Path
# Description: Simulates a single asset price path using Geometric Brownian Motion
# ensuring log-normal dynamics and positive prices.
# ==============================================================================

set.seed(123)

S0 <- 100
mu <- 0.05
sigma <- 0.2
T <- 1
n_steps <- 100
dt <- T / n_steps

S <- numeric(n_steps)
S[1] <- S0

Z <- rnorm(n_steps)

for (i in 2:n_steps) {
  S[i] <- S[i-1] * exp((mu - 0.5 * sigma^2) * dt +
                         sigma * sqrt(dt) * Z[i])
}

plot(S, type = "l", main = "Simulated GBM Path",
     xlab = "Time Step", ylab = "Price")


# ==============================================================================
# SECTION 6: Multiple GBM Price Paths
# Description: Generates multiple independent asset price paths using GBM 
# and visualizes them together on a single plot.
# ==============================================================================

set.seed(123)

n_paths <- 50
n_steps <- 100
S_paths <- matrix(0, nrow = n_steps, ncol = n_paths)

for (j in 1:n_paths) {
  Z <- rnorm(n_steps)
  S <- numeric(n_steps)
  S[1] <- 100
  
  for (i in 2:n_steps) {
    S[i] <- S[i-1] * exp((0.05 - 0.5 * 0.2^2) * (1/n_steps) +
                           0.2 * sqrt(1/n_steps) * Z[i])
  }
  
  S_paths[, j] <- S
}

matplot(S_paths, type = "l", lty = 1,
        main = "Multiple Simulated Price Paths",
        xlab = "Time Step", ylab = "Price")


# ==============================================================================
# SECTION 7: Payoff Evaluation
# Description: Computes the terminal payoff of a European call option for each 
# of the simulated paths from the previous step.
# ==============================================================================

K <- 100

S_T <- S_paths[n_steps, ]
payoffs <- pmax(S_T - K, 0)

head(payoffs)


# ==============================================================================
# SECTION 8: Estimation of Option Price
# Description: Discounts the expected terminal payoff back to present value 
# using the continuous risk-free rate to find the option price.
# ==============================================================================

r <- 0.05
T <- 1

price_estimate <- exp(-r * T) * mean(payoffs)
price_estimate


# ==============================================================================
# SECTION 9: Asian Option Pricing Functions
# Description: Defines functions to price Asian Call and Put options by 
# averaging the asset price over the simulated path.
# ==============================================================================

price_asian_call <- function(S0, K, r, sigma, T, steps, n_paths) {
  S <- simulate_paths(S0, r, sigma, T, steps, n_paths)
  
  avg_price <- rowMeans(S[, -1])
  payoff <- pmax(avg_price - K, 0)
  
  return(exp(-r * T) * mean(payoff))
}

price_asian_put <- function(S0, K, r, sigma, T, steps, n_paths) {
  S <- simulate_paths(S0, r, sigma, T, steps, n_paths)
  
  avg_price <- rowMeans(S[, -1])
  payoff <- pmax(K - avg_price, 0)
  
  return(exp(-r * T) * mean(payoff))
}


# ==============================================================================
# SECTION 10: Barrier Option Pricing Function
# Description: Defines a function to price barrier options (Up/Down & In/Out) 
# by checking if paths cross a defined barrier boundary.
# ==============================================================================

price_barrier_call <- function(S0, K, r, sigma, T, B,
                               barrier_type = "do",
                               steps, n_paths) {
  
  S <- simulate_paths(S0, r, sigma, T, steps, n_paths)
  
  breached <- apply(S, 1, function(path) {
    if (barrier_type == "do") return(any(path <= B))
    if (barrier_type == "di") return(any(path <= B))
    if (barrier_type == "uo") return(any(path >= B))
    if (barrier_type == "ui") return(any(path >= B))
  })
  
  ST <- S[, ncol(S)]
  payoff <- pmax(ST - K, 0)
  
  # Apply barrier logic
  if (barrier_type == "do" || barrier_type == "uo") {
    payoff[breached] <- 0
  } else if (barrier_type == "di" || barrier_type == "ui") {
    payoff[!breached] <- 0
  }
  
  return(exp(-r * T) * mean(payoff))
}


# ==============================================================================
# SECTION 11: Lookback Option Pricing Functions
# Description: Defines functions for pricing fixed and floating strike lookback 
# options by tracking the minimum or maximum price of a given path.
# ==============================================================================

price_lookback_call <- function(S0, K, r, sigma, T,
                                strike_type = "fixed",
                                steps, n_paths) {
  
  S <- simulate_paths(S0, r, sigma, T, steps, n_paths)
  
  if (strike_type == "fixed") {
    max_price <- apply(S, 1, max)
    payoff <- pmax(max_price - K, 0)
    
  } else if (strike_type == "floating") {
    min_price <- apply(S, 1, min)
    ST <- S[, ncol(S)]
    payoff <- pmax(ST - min_price, 0)
  }
  
  return(exp(-r * T) * mean(payoff))
}

price_lookback_put <- function(S0, K, r, sigma, T,
                               strike_type = "fixed",
                               steps, n_paths) {
  
  S <- simulate_paths(S0, r, sigma, T, steps, n_paths)
  
  if (strike_type == "fixed") {
    min_price <- apply(S, 1, min)
    payoff <- pmax(K - min_price, 0)
    
  } else if (strike_type == "floating") {
    max_price <- apply(S, 1, max)
    ST <- S[, ncol(S)]
    payoff <- pmax(max_price - ST, 0)
  }
  
  return(exp(-r * T) * mean(payoff))
}


# ==============================================================================
# SECTION 12: Generic Path Simulation Wrapper
# Description: The core engine utilized by the exotic pricing functions above 
# to cleanly generate matrix outputs of multiple GBM paths.
# ==============================================================================

simulate_paths <- function(S0, r, sigma, T, steps, n_paths) {
  dt <- T / steps
  S  <- matrix(0, nrow = n_paths, ncol = steps + 1)
  S[, 1] <- S0
  for (i in 2:(steps + 1)) {
    Z      <- rnorm(n_paths)
    S[, i] <- S[, i-1] * exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
  }
  return(S)
}


# ==============================================================================
# SECTION 13: Black-Scholes Baselines and Setup
# Description: Defines the true Black-Scholes analytical formulas for European 
# calls/puts to serve as benchmarks for the Monte Carlo estimations.
# ==============================================================================

set.seed(42)

S0      <- 100
K       <- 100
r       <- 0.05
sigma   <- 0.2
T       <- 1
steps   <- 252
n_paths <- 10000
B       <- 85

bs_call <- function(S0, K, r, sigma, T) {
  d1 <- (log(S0/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  S0*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
}

bs_put <- function(S0, K, r, sigma, T) {
  d1 <- (log(S0/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  K*exp(-r*T)*pnorm(-d2) - S0*pnorm(-d1)
}

bs_call_price <- bs_call(S0, K, r, sigma, T)
bs_put_price  <- bs_put(S0, K, r, sigma, T)


# ==============================================================================
# SECTION 14: Comprehensive Option Pricing Outputs
# Description: Executes the simulation over all defined option types and 
# compiles the results vs theoretical benchmarks into a single data frame.
# ==============================================================================

# European: direct simulate_paths (no Section 5 wrapper exists)
S       <- simulate_paths(S0, r, sigma, T, steps, n_paths)
ST      <- S[, steps + 1]
eu_call <- exp(-r*T) * mean(pmax(ST - K, 0))
eu_put  <- exp(-r*T) * mean(pmax(K - ST, 0))

# All exotic options via Section 5 functions
asian_call    <- price_asian_call(S0, K, r, sigma, T, steps, n_paths)
asian_put     <- price_asian_put(S0, K, r, sigma, T, steps, n_paths)

barrier_do    <- price_barrier_call(S0, K, r, sigma, T,
                                    B = B, barrier_type = "do",
                                    steps = steps, n_paths = n_paths)

lb_call_fixed <- price_lookback_call(S0, K, r, sigma, T,
                                     strike_type = "fixed",
                                     steps = steps, n_paths = n_paths)
lb_put_fixed  <- price_lookback_put(S0, K, r, sigma, T,
                                    strike_type = "fixed",
                                    steps = steps, n_paths = n_paths)


results <- data.frame(
  Option   = c("European Call", "European Put",
               "Asian Call", "Asian Put",
               "Barrier (Down-and-Out Call)",
               "Lookback Call (Fixed Strike)",
               "Lookback Put (Fixed Strike)"),
  MC_Price = round(c(eu_call, eu_put, asian_call, asian_put,
                     barrier_do, lb_call_fixed, lb_put_fixed), 4),
  BS_Price = round(c(bs_call_price, bs_put_price,
                     NA, NA, NA, NA, NA), 4)
)
print(results)


# ==============================================================================
# SECTION 15: Convergence Analysis by Number of Paths (N)
# Description: Examines standard error decay against theoretical rates as N grows, 
# plotting price convergence alongside the O(N^-1/2) error curve.
# ==============================================================================

set.seed(42)

path_sizes <- c(100, 500, 1000, 2000, 5000, 10000, 25000, 50000)
mc_eu      <- numeric(length(path_sizes))
mc_eu_se   <- numeric(length(path_sizes))

for (idx in seq_along(path_sizes)) {
  N             <- path_sizes[idx]
  S_conv        <- simulate_paths(S0, r, sigma, T, steps, N)
  pv_eu         <- exp(-r*T) * pmax(S_conv[, steps + 1] - K, 0)
  mc_eu[idx]    <- mean(pv_eu)
  mc_eu_se[idx] <- sd(pv_eu) / sqrt(N)
}

price_range <- range(c(mc_eu, bs_call_price))
price_pad   <- diff(price_range) * 0.2
se_ref      <- mc_eu_se[1] * sqrt(path_sizes[1]) / sqrt(path_sizes)
se_range    <- range(c(mc_eu_se, se_ref))

par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))

# Plot 1: Price convergence
plot(path_sizes, mc_eu,
     type = "b", log = "x", pch = 16, col = "black",
     xlab = "Number of Paths N (log scale)",
     ylab = "Estimated Call Price",
     main = "MC Price Convergence vs N\n(European Call)",
     ylim = c(price_range[1] - price_pad, price_range[2] + price_pad),
     xaxt = "n")
axis(1, at = path_sizes,
     labels = format(path_sizes, big.mark = ",", scientific = FALSE),
     las = 2, cex.axis = 0.7)
abline(h = seq(price_range[1] - price_pad,
               price_range[2] + price_pad, length.out = 8),
       col = "grey90", lty = 3, lwd = 0.8)
abline(h = bs_call_price, col = "red", lwd = 2, lty = 2)
lines(path_sizes, mc_eu, type = "b", pch = 16)
legend("topright",
       legend = c("MC Estimate", "BS Benchmark"),
       col = c("black", "red"), lty = c(1, 2),
       lwd = 2, pch = c(16, NA), bty = "n", cex = 0.85)

# Plot 2: SE decay (log-log)
plot(path_sizes, mc_eu_se,
     type = "b", log = "xy", pch = 16, col = "black",
     xlab = "Number of Paths N (log scale)",
     ylab = "Standard Error (log scale)",
     main = "Standard Error Decay vs N\n(log-log scale)",
     ylim = c(min(se_range) * 0.7, max(se_range) * 1.4),
     xaxt = "n")
axis(1, at = path_sizes,
     labels = format(path_sizes, big.mark = ",", scientific = FALSE),
     las = 2, cex.axis = 0.7)
lines(path_sizes, se_ref, col = "steelblue", lty = 2, lwd = 2)
legend("topright",
       legend = c("Empirical SE", expression(O(N^{-1/2}))),
       col = c("black", "steelblue"), lty = c(1, 2),
       lwd = 2, pch = c(16, NA), bty = "n", cex = 0.85)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))


# ==============================================================================
# SECTION 16: Confidence Interval Analysis
# Description: Constructs 95% CI data frames and plots absolute error and interval 
# width to visually confirm statistical bias and estimation stability.
# ==============================================================================

errors <- abs(mc_eu - bs_call_price)

ci_data <- data.frame(
  N       = path_sizes,
  Price   = round(mc_eu, 4),
  SE      = round(mc_eu_se, 5),
  CI_low  = round(mc_eu - 1.96*mc_eu_se, 4),
  CI_high = round(mc_eu + 1.96*mc_eu_se, 4),
  Width   = round(2*1.96*mc_eu_se, 5),
  Abs_Err = round(errors, 5)
)
print(ci_data)

# Anchor reference line on median to avoid near-zero instability
ref_anchor <- median(errors) * sqrt(median(path_sizes))
err_ref    <- ref_anchor / sqrt(path_sizes)
err_ylim   <- c(min(errors[errors > 0]) * 0.4,
                max(errors) * 2.5)
width_ylim <- c(0, max(ci_data$Width) * 1.2)

par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))

# Plot 1: Absolute error log-log
plot(path_sizes, errors,
     type = "b", log = "xy", pch = 16, col = "black",
     xlab = "Number of Paths N (log scale)",
     ylab = "|MC Price - BS Price| (log scale)",
     main = "Absolute Error vs N\n(log-log scale)",
     ylim = err_ylim,
     xaxt = "n")
axis(1, at = path_sizes,
     labels = format(path_sizes, big.mark = ",", scientific = FALSE),
     las = 2, cex.axis = 0.7)
lines(path_sizes, err_ref, col = "red", lty = 2, lwd = 2)
legend("topright",
       legend = c("Absolute Error", expression(O(N^{-1/2}))),
       col = c("black", "red"), lty = c(1, 2),
       lwd = 2, pch = c(16, NA), bty = "n", cex = 0.85)

# Plot 2: CI width
plot(path_sizes, ci_data$Width,
     type = "b", log = "x", pch = 16, col = "black",
     xlab = "Number of Paths N (log scale)",
     ylab = "95% Confidence Interval Width",
     main = "CI Width vs N",
     ylim = width_ylim,
     xaxt = "n")
axis(1, at = path_sizes,
     labels = format(path_sizes, big.mark = ",", scientific = FALSE),
     las = 2, cex.axis = 0.7)
abline(h = seq(0, max(ci_data$Width) * 1.2, length.out = 7),
       col = "grey90", lty = 3, lwd = 0.8)
lines(path_sizes, ci_data$Width, type = "b", pch = 16)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))


# ==============================================================================
# SECTION 17: Convergence by Time Discretization (Steps)
# Description: Evaluates exotic option price sensitivity to path coarseness by 
# varying the number of internal time steps.
# ==============================================================================

set.seed(42)

step_sizes      <- c(10, 25, 50, 100, 252, 500)
barrier_prices  <- numeric(length(step_sizes))
lookback_prices <- numeric(length(step_sizes))
asian_prices    <- numeric(length(step_sizes))

for (idx in seq_along(step_sizes)) {
  m <- step_sizes[idx]
  barrier_prices[idx]  <- price_barrier_call(S0, K, r, sigma, T,
                                             B = B, barrier_type = "do",
                                             steps = m, n_paths = n_paths)
  lookback_prices[idx] <- price_lookback_call(S0, K, r, sigma, T,
                                              strike_type = "fixed",
                                              steps = m, n_paths = n_paths)
  asian_prices[idx]    <- price_asian_call(S0, K, r, sigma, T,
                                           steps = m, n_paths = n_paths)
}

# Padded ylim helper
ylim_pad <- function(x, frac = 0.15) {
  rng <- range(x)
  pad <- diff(rng) * frac
  if (pad < 1e-6) pad <- abs(rng[1]) * 0.05
  c(rng[1] - pad, rng[2] + pad)
}

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

plot(step_sizes, barrier_prices,
     type = "b", pch = 16, col = "black",
     xlab = "Number of Time Steps",
     ylab = "Estimated Price",
     main = "Barrier (Down-and-Out Call)\nvs Time Steps",
     ylim = ylim_pad(barrier_prices),
     xaxt = "n")
axis(1, at = step_sizes)
abline(h = pretty(barrier_prices, n = 6),
       col = "grey90", lty = 3, lwd = 0.8)
lines(step_sizes, barrier_prices, type = "b", pch = 16)

plot(step_sizes, lookback_prices,
     type = "b", pch = 16, col = "black",
     xlab = "Number of Time Steps",
     ylab = "Estimated Price",
     main = "Lookback Call (Fixed Strike)\nvs Time Steps",
     ylim = ylim_pad(lookback_prices),
     xaxt = "n")
axis(1, at = step_sizes)
abline(h = pretty(lookback_prices, n = 6),
       col = "grey90", lty = 3, lwd = 0.8)
lines(step_sizes, lookback_prices, type = "b", pch = 16)

plot(step_sizes, asian_prices,
     type = "b", pch = 16, col = "black",
     xlab = "Number of Time Steps",
     ylab = "Estimated Price",
     main = "Asian Call\nvs Time Steps",
     ylim = ylim_pad(asian_prices),
     xaxt = "n")
axis(1, at = step_sizes)
abline(h = pretty(asian_prices, n = 6),
       col = "grey90", lty = 3, lwd = 0.8)
lines(step_sizes, asian_prices, type = "b", pch = 16)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))