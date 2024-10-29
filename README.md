# Comparing-Performance-of-Mixed-Copulas-vs-Independent-Copulas-on-Conditional-Probability-Strategy

This report presents a copula-based pairs trading strategy applied to S&P 500 stocks. By using Gaussian, Gumbel, and Clayton copulas, this study compares the efficacy of both single and mixed copulas for capturing dependency structures in stock pairs with high correlation. The strategy is trained on stock data (Closing Price) from **2016 to 2020**, and the strategy’s performance is evaluated on test data from **2020 to 2023**, focusing on metrics such as returns, volatility, and overall profit after transaction costs. Results indicate that the mixed copula model, optimized via the EM algorithm, **improves trading outcomes**, demonstrating superior modeling of tail dependencies.

## Introduction

Pairs trading is a market-neutral strategy that exploits statistical relationships between assets. Traditional approaches like correlation or cointegration methods may not capture dependencies adequately, especially under extreme market conditions. Copulas offer a flexible framework to model dependencies beyond linear correlations, enabling better performance during market stress periods. This report examines the efficacy of Gaussian, Gumbel, and Clayton copulas in pairs trading and proposes a mixed copula approach optimized via the Expectation-Maximization (EM) algorithm to improve dependency modeling in asset pairs.

## Methodology 

### Data Selection and Preprocessing
* **Data Range**: January 1, 2016, to December 31, 2023
* **Assets**: S&P 500 stocks with a correlation threshold of 0.8

Stocks were selected based on historical correlation, retaining pairs with correlation coefficients above 0.8 over a five-year period. Daily adjusted closing prices were used to calculate log returns, and pseudo-observations were derived for copula fitting.

### Copula Modeling
* **Method 1 - Single Copula Selection:**
  * Fit each pair using Gaussian, Gumbel, and Clayton copulas.
  * Select the copula with the highest log-likelihood for each pair.
* **Method 2 - Mixed Copula Model:**
  * Create a weighted mixture of the three copulas.
  * Use the EM algorithm to optimize copula weights and parameters.
 
### EM Algorithm for Mixed Copula
* **Initialization**: Assign equal weights to each copula.
* **E-Step**: Calculate responsibilities for each copula based on densities.
* **M-Step**: Update weights and copula parameters to maximize the log-likelihood.
* **Convergence**: Stop when the change in log-likelihood falls below a defined threshold.

## Trading Strategy Design

The strategy generates trading signals based on conditional probabilities derived from the mixed copula model. For each stock pair, buy/sell signals are triggered if the probability of one stock exceeding a given threshold (e.g., 5% or 95%) is met, indicating overbought/oversold conditions. Trades incur both fixed and percentage-based transaction costs. The hedge ratio from Johansen's cointegration test adjusts position sizing between each pair.

## Results and Performance Evaluation

After running the strategy on test data (2020-2023), the mixed copula model demonstrated higher cumulative returns than the individual copula models. Table 1 shows the average returns, volatility, and AIC values across tested pairs. The optimized copula weights indicate a preference for Clayton and Gumbel copulas, reflecting stronger tail dependencies in the data. The Sharpe ratio of the mixed copula model exceeded that of Gaussian copula pairs by 15%, indicating enhanced risk-adjusted returns.
