library(quantmod)
library(copula)
library(urca)

find_correlation <- function(stock1, stock2, start_date, end_date) {
  adjusted_col1 <- paste0(stock1, ".Close")
  adjusted_col2 <- paste0(stock2, ".Close")
  
  # Fetching the stock data with error handling
  df1 <- tryCatch({
    as.data.frame(getSymbols(stock1, src='yahoo', auto.assign=FALSE, from=start_date, to=end_date))[[adjusted_col1]]
  }, error = function(e) {
    warning(paste("Failed to fetch data for", stock1))
    return(NULL)
  })
  
  df2 <- tryCatch({
    as.data.frame(getSymbols(stock2, src='yahoo', auto.assign=FALSE, from=start_date, to=end_date))[[adjusted_col2]]
  }, error = function(e) {
    warning(paste("Failed to fetch data for", stock2))
    return(NULL)
  })
  
  # Check if either df1 or df2 is NULL
  if (is.null(df1) || is.null(df2)) {
    return(NA)  # Return NA if data fetching failed
  }
  
  # Ensure the lengths match before calculating correlation
  min_length <- min(length(df1), length(df2))
  correlation <- cor(df1[1:min_length], df2[1:min_length], use="complete.obs")
  
  return(correlation)
}

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

copula_conditional <- function(stocks, train_start, train_end, test_start, test_end, p1, p2, fixed_cost, percentage_cost){
  # training data
  adjusted_col1 <- paste0(stocks[1], ".Close")
  adjusted_col2 <- paste0(stocks[2], ".Close")
  
  df_train_1 <- as.data.frame(getSymbols(stocks[1], src='yahoo', auto.assign=FALSE, from=train_start, to=train_end))[adjusted_col1]
  df_train_2 <- as.data.frame(getSymbols(stocks[2], src='yahoo', auto.assign=FALSE, from=train_start, to=train_end))[adjusted_col2]
  
  df_train <- data.frame(s1 = df_train_1[[adjusted_col1]], s2 = df_train_2[[adjusted_col2]])  
  
  colnames(df_train) <- c(stocks[1], stocks[2])
  
  # Calculate log returns
  #returns <- diff(log(df_train))[-1, ]
  #returns <- as.data.frame(returns)  # Convert to data frame for compatibility
  
  convert_to_pobs <- function(data) {
    # Convert data to matrix
    data_matrix <- as.matrix(data)
    
    # Convert data to pseudo-observations
    pobs_data <- pobs(data_matrix)
    
    # Check if pseudo-observations contain NA values
    if (any(is.na(pobs_data))) {
      stop("Pseudo-observations contain NA values.")
    }
    
    return(pobs_data)
  }
  u <- convert_to_pobs(df_train)
  
  weights <- em(stocks, train_start, train_end)$optimized_weights
  
  # Mixed Conditional Probability Function
  combined_conditional_probability <- function(x, y, variable, weights) {
    if(variable ==1){
      conditional_gaussian <- pCopula(cbind(x, y), gaussian_fit$copula) / pCopula(c(1, mat_test[i, 2]), gaussian_fit$copula)
      conditional_gumbel <- pCopula(cbind(x, y), gumbel_fit$copula) / pCopula(c(1, mat_test[i, 2]), gumbel_fit$copula)
      conditional_clayton <- pCopula(cbind(x, y), clayton_fit$copula) / pCopula(c(1, mat_test[i, 2]), clayton_fit$copula)    
    }
    else{
      conditional_gaussian <- pCopula(cbind(x, y), gaussian_fit$copula) / pCopula(c(mat_test[i, 1], 1), gaussian_fit$copula)
      conditional_gumbel <- pCopula(cbind(x, y), gumbel_fit$copula) / pCopula(c(mat_test[i, 1], 1), gumbel_fit$copula)
      conditional_clayton <- pCopula(cbind(x, y), clayton_fit$copula) / pCopula(c(mat_test[i, 1], 1), clayton_fit$copula)
    }
    return(weights[1] * conditional_gaussian + weights[2] * conditional_gumbel + weights[3] * conditional_clayton)
  }
  
  
  # testing data
  
  df_test_1 = as.data.frame(getSymbols(stocks[1], src='yahoo', auto.assign=FALSE, from=test_start, to=test_end))[adjusted_col1]
  df_test_2 = as.data.frame(getSymbols(stocks[2], src='yahoo', auto.assign=FALSE, from=test_start, to=test_end))[adjusted_col2]
  
  df_test <- data.frame(s1 = df_test_1[[adjusted_col1]], s2 = df_test_2[[adjusted_col2]])
  
  # Johansen cointegration test on training data
  jtest_train <- ca.jo(df_train, type="trace", K=2, ecdet="none", spec="longrun")
  
  # Extract cointegration vectors and calculate hedge ratio
  cointegration_vectors_train <- jtest_train@V
  hedge_ratio <- cointegration_vectors_train[1, 1] / cointegration_vectors_train[2, 1]
  
  # Calculate empirical CDF values for testing data
  u_s1_test <- ecdf(df_test$s1)(df_test$s1)
  u_s2_test <- ecdf(df_test$s2)(df_test$s2)
  
  # Prepare matrix for copula fitting on testing data
  mat_test <- cbind(u_s1_test, u_s2_test)
  
  # Initialize variables for trading strategy
  initial_capital <- 100  # Starting capital
  money <- initial_capital
  position <- 0
  equity <- c()  # Vector to store returns for each trading day
  
  # Trading logic on testing data
  for (i in 1:nrow(df_test)) {
    
    p_copula_s1 <- combined_conditional_probability(mat_test[i,1], mat_test[i,2], 1, weights)
    p_copula_s2 <- combined_conditional_probability(mat_test[i,1], mat_test[i,2], 2, weights)
    
    
    trade_cost <- 0
    
    if (p_copula_s1 <= p1 | p_copula_s2 >= p2 & position == 0) {
      money <- money - (hedge_ratio * df_test$s1[i]) + (df_test$s2[i])
      trade_value <- abs(hedge_ratio * df_test[,1][i]) + abs(df_test[,2][i])
      trade_cost <- fixed_cost + (percentage_cost * trade_value)
      money <- money - trade_cost
      if(money!=0){
        equity <- c(equity, money)
      }
      position <- 1
      
    } else if (p_copula_s1 >= p2 & p_copula_s2 <= p1 & position == 1) {
      money <- money + (hedge_ratio * df_test$s1[i]) - (df_test$s2[i])
      trade_value <- abs(hedge_ratio * df_test[,1][i]) + abs(df_test[,2][i])
      trade_cost <- fixed_cost + (percentage_cost * trade_value)
      money <- money - trade_cost
      if(money!=0){
        equity <- c(equity, money)
      }
      position <- 0
    }
    
  }
  
  # Print final result and Sharpe Ratio
  eq <- c()
  for (i in 1:length(equity)){
    if(equity[i]!=0){
      eq <- c(eq, equity[i])
    }
  }
  
  ret <- c()
  for (i in 2:length(eq)){
    ret <- c(ret, ((eq[i]-eq[i-1])/eq[i-1]))
  }
  
  volatility <- sd(ret)
  
  return (list(money = money, volatility=volatility))
}

stocks <- c(
  "AAPL", "MSFT", "AMZN", "GOOGL", "FB", "TSLA", "JNJ", "V", "NVDA",
  "PYPL", "UNH", "HD", "PG", "DIS", "MA", "NFLX", "VZ", "ADBE", "CMCSA",
  "PFE", "INTC", "KO", "MRK", "PEP", "ABBV", "T", "CSCO", "XOM", "WMT",
  "NKE", "CVX", "LLY", "MCD", "QCOM", "BMY", "ACN", "AVGO", "HON", "COST",
  "IBM", "MDT", "TXN", "AMGN", "DHR", "NEE", "LIN", "UNP"
)

start_date <- "2016-01-01"
end_date <- "2020-01-01"
pairs <- data.frame(stock1 = character(), stock2 = character(), correlation = numeric())

# Generate unique stock pairs
combinations <- combn(stocks, 2, simplify = FALSE)

for (pair in combinations) {
  stock1 <- pair[1]
  stock2 <- pair[2]

  correlation <- find_correlation(stock1, stock2, start_date, end_date)

  if (!is.na(correlation) && correlation >= 0.8) {
    pairs <- rbind(pairs, data.frame(stock1 = stock1, stock2 = stock2, correlation = round(correlation, 4)))
  }
}

selected_pairs <- matrix(, nrow=nrow(pairs), ncol=2)

for (i in 1:nrow(pairs)){
  selected_pairs[i, 1]=pairs$stock1[i]
  selected_pairs[i, 2]=pairs$stock2[i]
}

train_start <- '2015-01-01'
train_end <-'2020-01-01'
test_start <- '2021-01-01'
test_end <- '2023-12-31'
p1 <- 0.05
p2 <- 0.95
fixed_cost <- 10
percentage_cost <- 0.001

final_result_list <- list()  # Initialize a list to store results
vis <- c()

for (i in 1:nrow(selected_pairs)){
  final <- copula_conditional(c(selected_pairs[i, 1], selected_pairs[i, 2]), train_start, train_end, test_start, test_end, p1, p2, fixed_cost, percentage_cost)
  # Store the results in a temporary data frame
  temp_result <- data.frame(
    Stock1 = selected_pairs[i, 1],
    Stock2 = selected_pairs[i, 2],
    Returns = final$money,
    Volatility = final$volatility
  )
  
  # Append the temporary result to the list
  final_result_list[[length(final_result_list) + 1]] <- temp_result
}

# Combine all results into a single data frame
final_result_df <- do.call(rbind, final_result_list)

# Export as csv
write.csv(final_result_df, file = "Copula_results", row.names = FALSE)

# Print the final results DataFrame
print(final_result_df)
