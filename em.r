library(quantmod)

em <- function(tickers, start_date, end_date) {
  # Step 1: Fetch stock data from Yahoo Finance
  getSymbols(tickers, src = "yahoo", from = start_date, to = end_date)
  
  # Merge and extract closing prices
  stock_data <- merge(Cl(get(tickers[1])), Cl(get(tickers[2])))
  colnames(stock_data) <- tickers
  
  # Step 2: Calculate log returns
  calculate_log_returns <- function(data) {
    log_returns <- diff(log(data))
    return(as.data.frame(log_returns[-1, ]))  # Remove the first NA row
  }
  
  log_returns <- calculate_log_returns(stock_data)
  
  # Transform returns to uniform margins using pobs (pseudo-observations)
  u <- pobs(log_returns)
  
  # Step 3: Initialize the copulas
  gaussian_cop <- normalCopula(dim = 2)
  gumbel_cop <- gumbelCopula(dim = 2)
  clayton_cop <- claytonCopula(dim = 2)
  
  # Initial weights for each copula
  initial_weights <- c(1/3, 1/3, 1/3)
  
  # Initial parameter estimates for each copula based on the uniform data
  initial_params <- list(
    gaussian = fitCopula(gaussian_cop, u, method = "ml")@estimate,
    gumbel = fitCopula(gumbel_cop, u, method = "ml")@estimate,
    clayton = fitCopula(clayton_cop, u, method = "ml")@estimate
  )
  
  # Step 4: Function to calculate the density for each copula
  copula_density <- function(u, copula, params) {
    copula <- setTheta(copula, params)  # Set the copula parameters
    return(dCopula(u, copula))
  }
  
  # Step 5: E-step - Calculate responsibilities (gamma values)
  calc_responsibilities <- function(u, weights, params) {
    n <- nrow(u)
    K <- length(weights)  # Number of copulas
    gamma <- matrix(0, n, K)
    
    # Densities for each copula
    densities <- list(
      gaussian = copula_density(u, gaussian_cop, params$gaussian),
      gumbel = copula_density(u, gumbel_cop, params$gumbel),
      clayton = copula_density(u, clayton_cop, params$clayton)
    )
    
    # Calculate responsibilities
    for (k in 1:K) {
      gamma[, k] <- weights[k] * densities[[k]]
    }
    
    # Normalize to make them probabilities
    gamma <- gamma / rowSums(gamma)
    return(gamma)
  }
  
  # Step 6: M-step - Update weights and parameters
  update_parameters <- function(u, gamma) {
    new_weights <- colMeans(gamma)
    
    new_params <- list()
    new_params$gaussian <- fitCopula(gaussian_cop, u, weights = gamma[, 1], method = "ml")@estimate
    new_params$gumbel <- fitCopula(gumbel_cop, u, weights = gamma[, 2], method = "ml")@estimate
    new_params$clayton <- fitCopula(clayton_cop, u, weights = gamma[, 3], method = "ml")@estimate
    
    return(list(weights = new_weights, params = new_params))
  }
  
  # Step 7: Function to calculate the log-likelihood of the mixed copula model
  calc_log_likelihood <- function(u, weights, params) {
    densities <- list(
      gaussian = copula_density(u, gaussian_cop, params$gaussian),
      gumbel = copula_density(u, gumbel_cop, params$gumbel),
      clayton = copula_density(u, clayton_cop, params$clayton)
    )
    
    # Mixed density based on weights
    mixed_density <- weights[1] * densities$gaussian + 
      weights[2] * densities$gumbel + 
      weights[3] * densities$clayton
    
    # Log-likelihood
    log_likelihood <- sum(log(mixed_density))
    return(log_likelihood)
  }
  
  # Step 8: EM algorithm implementation
  em_algorithm <- function(u, initial_weights, initial_params, tol = 1e-6, max_iter = 100) {
    weights <- initial_weights
    params <- initial_params
    log_likelihoods <- c()
    
    for (iter in 1:max_iter) {
      # E-step
      gamma <- calc_responsibilities(u, weights, params)
      
      # M-step
      updates <- update_parameters(u, gamma)
      weights <- updates$weights
      params <- updates$params
      
      # Calculate log-likelihood
      log_likelihood <- calc_log_likelihood(u, weights, params)
      log_likelihoods <- c(log_likelihoods, log_likelihood)
      
      # Check for convergence
      if (iter > 1 && abs(log_likelihoods[iter] - log_likelihoods[iter - 1]) < tol) {
        cat("Convergence reached at iteration:", iter, "\n")
        break
      }
    }
    
    return(list(weights = weights, params = params, log_likelihood = log_likelihoods))
  }
  
  # Run the EM algorithm
  result <- em_algorithm(u, initial_weights, initial_params)
  
  # Log-likelihood for independent copula (Gaussian)
  independent_log_likelihood <- sum(dCopula(u, normalCopula(dim = 2), log = TRUE))
  
  # Calculate final AIC using the last log-likelihood and optimized weights
  final_log_likelihood <- calc_log_likelihood(u, result$weights, result$params)
  num_params <- 3  # Three weights to estimate
  final_aic <- -2 * final_log_likelihood + 2 * num_params
  
  # Return results
  return(list(
    optimized_weights = result$weights,
    final_params = result$params,
    final_log_likelihood = final_log_likelihood,
    independent_log_likelihood = independent_log_likelihood,
    final_aic = final_aic
  ))
}
