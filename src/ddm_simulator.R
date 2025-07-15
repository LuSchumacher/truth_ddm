sample_ddm <- function(v = 1, a = 1, ndt = 0.3, bias = 0.5) {
  # Drift Diffusion Model (DDM) sampler
  # v     = drift rate
  # a     = boundary separation
  # ndt   = non-decision time (in seconds)
  # bias  = starting point bias (0.5 = unbiased)
  
  max_iter <- 1e4         # max number of time steps
  dt <- 0.001             # time step size
  noise_sd <- 1           # noise standard deviation
  sqrt_dt <- sqrt(dt)     # precompute for efficiency
  
  # Initialize path and generate noise
  z <- a * bias           # starting point
  path <- z
  noise <- rnorm(max_iter, mean = 0, sd = noise_sd)
  
  # Accumulate evidence
  iter <- 1
  while (path > 0 && path < a && iter < max_iter) {
    path <- path + v * dt + sqrt_dt * noise[iter]
    iter <- iter + 1
  }
  
  # Determine decision (0 = lower boundary, 1 = upper boundary)
  decision <- as.numeric(path >= a)
  rt <- ndt + iter * dt
  
  return(c(decision, rt))
}