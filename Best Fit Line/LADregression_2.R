ladfit <- function(x, y) {
  if (length(x) != length(y)) stop("Lengths and x and y do not match")
  skip <- is.na(x) | is.na(y)
  x <- x[!skip]
  y <- y[!skip]
  lambda <- function(par) {
    a <- par[1]
    b <- par[2]
    sum(abs(y - a - b * x))
  }
  opt <- optim(par = c(0, 0), fn = lambda)
  opt$par
}

x <- rnorm(50)
y <- rnorm(50) + x
ladfit(x,y)
plot(x,y, main = "LAD regression", pch = 16)
abline(a = ladfit(x,y)[1], b = ladfit(x,y)[2], lwd = 2 , col = 'red')


# My approach....

theta_values <- seq(0, 180, by = 1)          # creating values of theta
theta_a_hat_d_a <- data.frame(theta = theta_values, a_hat = NA, d_a = NA)   # creating the data frame


for(t in theta_values){
  
  theta_rad <- t*(pi/180)
  x_new <- x * cos(theta_rad) - y * sin(theta_rad)
  #y_new <- x * sin(theta_rad) + y * cos(theta_rad)     # this will create a cloud of data for one t
  
  incr = 0.1
  a = sort(x_new)[1]
  
  vec_total_devia <- c()
  vec_a <- c()
      
  while(a <= sort(x_new,decreasing = TRUE)[1]){
    
    total_devia <- sum(abs(a - x_new))
    vec_total_devia <- c(total_devia, vec_total_devia)
    vec_a <- c(a, vec_a)
   
    a = a + incr    
  }
  
  theta_a_hat_d_a[theta_a_hat_d_a$theta == t, "d_a"] <- sort(vec_total_devia)[1]
  theta_a_hat_d_a[theta_a_hat_d_a$theta == t, "a_hat"] <- vec_a[vec_total_devia == sort(vec_total_devia)[1]]
    # error at 50     
}


best_d_a <- sort(theta_a_hat_d_a$d_a)[1]
best_a_hat <- theta_a_hat_d_a$a_hat[theta_a_hat_d_a$d_a == sort(theta_a_hat_d_a$d_a)[1]]
best_theta <- theta_a_hat_d_a$theta[theta_a_hat_d_a$d_a == sort(theta_a_hat_d_a$d_a)[1]] 


slope <- tan((90 - best_theta) * (pi / 180))
y_intercept <- best_a_hat * sin(best_theta * (pi / 180)) - slope * best_a_hat * cos(best_theta * (pi / 180))
abline(y_intercept, slope)


