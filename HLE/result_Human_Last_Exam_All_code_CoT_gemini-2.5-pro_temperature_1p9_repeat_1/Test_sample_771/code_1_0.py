import numpy as np

# --- Parameters ---
initial_price = 100.0
annual_growth_rate = 0.10
# The problem states weekly moves of 20%, which we model as weekly volatility
weekly_volatility = 0.20
years_to_hold = 30
liquidation_years = 1
total_years = years_to_hold + liquidation_years
weeks_per_year = 52
num_simulations = 20000 # Number of simulated paths for Basket B

# --- Time setup ---
dt = 1.0 / weeks_per_year # time step in years
total_weeks = total_years * weeks_per_year
liquidation_start_week = years_to_hold * weeks_per_year

# --- Strategy A: Deterministic Growth ---
# To maximize value, sell at the end of the 31st year.
final_value_A = initial_price * (1.0 + annual_growth_rate)**total_years

# --- Strategy B: Monte Carlo Simulation for Volatile Growth ---

# We need to find the weekly drift (mu_w) that results in an *expected*
# annual growth of 10%.
# Expected weekly growth factor E[P(t)/P(t-1)] must equal the deterministic weekly growth factor.
deterministic_weekly_growth = (1.0 + annual_growth_rate)**dt
# The expected return for a log-normal process is exp(mu + sigma^2 / 2)
# E[P(t)/P(t-1)] = exp(mu_w + weekly_volatility**2 / 2)
# Solving for mu_w:
# mu_w + weekly_volatility**2 / 2 = log(deterministic_weekly_growth)
mu_w = np.log(deterministic_weekly_growth) - (weekly_volatility**2 / 2.0)

realized_values_B = []

print("Running Monte Carlo simulation for Basket B...")

for i in range(num_simulations):
    # Initialize the price path for this simulation
    prices_B = [initial_price]
    current_price = initial_price
    
    # Generate the price path for 31 years
    for _ in range(total_weeks):
        # Generate a random shock from a standard normal distribution
        z = np.random.normal(0, 1)
        # Calculate the return for this week and update the price
        # P(t) = P(t-1) * exp(mu_w + weekly_volatility * Z)
        current_price *= np.exp(mu_w + weekly_volatility * z)
        prices_B.append(current_price)

    # Liquidation strategy for B: find the maximum price in the 31st year
    # The 31st year corresponds to the last 52 weeks of the simulation.
    liquidation_period_prices = prices_B[liquidation_start_week:]
    max_value_in_liquidation = np.max(liquidation_period_prices)
    realized_values_B.append(max_value_in_liquidation)

# Calculate the expected (average) value from the simulations
expected_value_B = np.mean(realized_values_B)

print("\n--- Simulation Results ---")
print(f"Final Value of Strategy A (Deterministic):")
# Here, we format the equation and the result as requested.
print(f"{initial_price:.2f} * (1 + {annual_growth_rate})**{total_years} = {final_value_A:.2f}")

print("\nExpected Value of Strategy B (Volatile with optimal liquidation):")
print(f"Average of {num_simulations} simulations = {expected_value_B:.2f}")

print("\n--- Conclusion ---")
if expected_value_B > final_value_A:
    print("The strategy of choosing basket B has a higher expected value than that of choosing basket A.")
else:
    print("The strategy of choosing basket A has a higher expected value than that of choosing basket B.")
