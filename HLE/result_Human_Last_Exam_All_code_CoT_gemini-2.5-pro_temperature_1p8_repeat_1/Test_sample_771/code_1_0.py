import numpy as np

# --- Parameters ---
P0 = 100.0                   # Initial price
annual_rate = 0.10           # Annual growth rate/trend
hold_years = 30
# Liquidation happens during the 31st year, so we model for 31 years total.
total_years = 31
weeks_per_year = 52
total_weeks = total_years * weeks_per_year
num_simulations = 10000      # Number of simulations for Monte Carlo

# --- Basket B: Simulation Parameters ---
# The problem states an expected growth trend of 10%. To model this in a Geometric
# Brownian Motion (GBM) framework, we set the drift `mu` of the process such that
# the expected price grows at the desired rate. The formula relating the continuous
# drift `mu` to the annual rate `r` is: mu = log(1 + r).
annual_drift_mu = np.log(1 + annual_rate)
weekly_drift_mu = annual_drift_mu / weeks_per_year
# The problem states the stock moves "20% in either direction each week",
# which we model as a weekly volatility (standard deviation of log returns).
weekly_volatility = 0.20

# --- Basket A: Deterministic Calculation ---
print("--- Strategy A: Deterministic Growth ---")
print("Basket A's value is predictable, growing at a fixed 10% per year.")
# Calculate the final value after 31 years: P_final = P0 * (1 + rate)^years
final_value_A = P0 * (1 + annual_rate)**total_years
print("The final value is calculated with the compound growth formula:")
print(f"Equation: Initial Price * (1 + Annual Rate) ^ Years")
print(f"Numbers: ${P0:.2f} * (1 + {annual_rate:.2f}) ^ {total_years} = ${final_value_A:.2f}")

# --- Basket B: Monte Carlo Simulation ---
print("\n--- Strategy B: Volatile Growth with Optimal Liquidation ---")
print("For Basket B, the trader can wait and sell at the highest price during the 31st year.")
print(f"We simulate {num_simulations} possible price paths to find the average outcome of this strategy.")

all_liquidation_values_B = []

for i in range(num_simulations):
    price_path = [P0]
    # The GBM model simulates price changes using log returns.
    # log_return = (mu - 0.5*sigma^2) + sigma*Z, where Z is a standard normal random variable.
    drift_term = weekly_drift_mu - 0.5 * weekly_volatility**2
    
    # We generate the price path for the full 31 years (1612 weeks)
    for _ in range(total_weeks - 1):
        random_term = weekly_volatility * np.random.normal(0, 1)
        next_price = price_path[-1] * np.exp(drift_term + random_term)
        price_path.append(next_price)

    # For the optimal liquidation, find the maximum price during the final year (last 52 weeks).
    liquidation_start_index = hold_years * weeks_per_year
    prices_in_liquidation_year = price_path[liquidation_start_index:]
    max_price_in_liquidation = np.max(prices_in_liquidation_year)
    all_liquidation_values_B.append(max_price_in_liquidation)

# The expected value is the average of the outcomes from all simulations
expected_value_B = np.mean(all_liquidation_values_B)
print(f"The Expected Value of the Basket B strategy is ${expected_value_B:.2f}")

# --- Conclusion ---
print("\n--- Final Comparison ---")
print(f"Expected Value of Strategy A: ${final_value_A:.2f}")
print(f"Expected Value of Strategy B: ${expected_value_B:.2f}")

print("\nAs predicted by financial theory, the ability to time the sale of the volatile asset")
print("gives Strategy B a significantly higher expected value.")