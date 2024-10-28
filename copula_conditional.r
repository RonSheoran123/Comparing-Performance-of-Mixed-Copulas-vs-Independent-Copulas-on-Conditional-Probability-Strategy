library(quantmod)
library(copula)
library(urca)

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
