library(quantmod)

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
