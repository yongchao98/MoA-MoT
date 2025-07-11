import numpy as np

def solve_trader_problem():
    """
    This script simulates and compares the final values of two investment strategies
    to determine which has a higher expected value.
    """

    # --- Parameters ---
    # We normalize the starting price at the beginning of the 31st year to 100.
    # The core of the problem lies in the value generated during the liquidation year.
    p_start_liquidation_year = 100.0

    # Financial parameters from the problem description
    annual_growth_rate = 0.10
    weekly_volatility = 0.20  # Given as 20% movement each week
    weeks_in_year = 52

    # Simulation parameters
    num_simulations = 20000

    # --- Strategy A: Deterministic Growth ---
    # With no volatility, the way to maximize value is to sell at the end of the year.
    # The value increases by the annual growth rate over the liquidation year.
    final_value_A = p_start_liquidation_year * (1 + annual_growth_rate)

    # --- Strategy B: Volatile Growth with Optimal Liquidation ---
    # First, calculate the expected weekly growth rate from the annual rate.
    # (1 + weekly_rate)^52 = 1 + annual_rate
    weekly_growth_rate = (1 + annual_growth_rate)**(1 / weeks_in_year) - 1

    # This list will store the maximum price achieved in each simulation run.
    liquidation_values_B = []

    # Use a fixed seed for reproducible results.
    np.random.seed(42)

    print("Running Monte Carlo simulation for Strategy B...")
    for i in range(num_simulations):
        current_price = p_start_liquidation_year
        weekly_prices_in_run = []
        for _ in range(weeks_in_year):
            # Model the weekly price change using a multiplicative random walk.
            # Price_t = Price_{t-1} * (1 + drift + random_shock)
            random_shock = np.random.normal(loc=0, scale=weekly_volatility)
            current_price *= (1 + weekly_growth_rate + random_shock)
            # Stocks cannot have a negative price.
            if current_price < 0:
                current_price = 0
            weekly_prices_in_run.append(current_price)

        # The rational trader sells at the highest price during the year.
        max_price_in_year = max(weekly_prices_in_run)
        liquidation_values_B.append(max_price_in_year)

    # The expected value of Strategy B is the average of all the recorded maximums.
    expected_value_B = np.mean(liquidation_values_B)

    # --- Output the Results and Conclusion ---
    print("\n--- Comparison of Investment Strategies ---")
    print("The simulation compares the certain outcome of Strategy A against the expected (average) outcome of Strategy B.")
    
    # Print the final equation for Strategy A
    print("\nFinal Value of Strategy A:")
    print(f"{p_start_liquidation_year:.2f} * (1 + {annual_growth_rate:.2f}) = {final_value_A:.2f}")

    # Print the result for Strategy B
    print("\nExpected Final Value of Strategy B:")
    print(f"Average of {num_simulations:,} simulated maximums ≈ {expected_value_B:.2f}")

    print("\n--- Conclusion ---")
    if expected_value_B > final_value_A:
        print("The expected value of Strategy B is higher than the final value of Strategy A.")
        print("This is because the high volatility of Basket B creates a valuable timing option during the year-long liquidation period.")
    else:
        print("The simulation did not show a higher expected value for Strategy B.")

solve_trader_problem()