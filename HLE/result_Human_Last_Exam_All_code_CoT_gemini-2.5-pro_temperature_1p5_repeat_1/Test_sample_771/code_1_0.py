import math

def analyze_baskets():
    """
    Analyzes and compares the expected values of two investment strategies.
    """
    initial_investment = 1.0  # Use a normalized initial value of 1
    growth_rate = 0.10        # 10% growth rate/trend for both
    liquidation_year = 31     # Liquidation occurs over the 31st year

    print("--- Investment Basket Analysis ---")
    print(f"Comparing two baskets over a period of {liquidation_year} years, both with a 10% growth parameter.\n")

    # --- Basket A Calculation ---
    # Basket A has a guaranteed 10% compound growth. Its value is certain.
    # To maximize value, it's sold at the end of the 31st year.
    # The final equation for its value is: P0 * (1 + r)^t
    
    value_a = initial_investment * (1 + growth_rate)**liquidation_year

    print("Basket A (Compound Growth):")
    print("The final value is deterministic due to compounding growth with no volatility.")
    print("Final Value = Initial Investment * (1 + Growth Rate) ^ Years")
    print(f"Final Value = {initial_investment:.1f} * (1 + {growth_rate}) ^ {liquidation_year}")
    print(f"Final Value = {value_a:.4f}\n")


    # --- Basket B Calculation ---
    # Basket B is described as an ARIMA model with drift, implying a linear growth in its expectation.
    # E[Value] = P0 * (1 + r * t).
    # The optimal liquidation strategy gives a value greater than the expected value at any fixed point.
    # We calculate the expected value at the end of the 31st year as a lower bound.
    
    expected_value_b_lower_bound = initial_investment * (1 + growth_rate * liquidation_year)

    print("Basket B (Linear Expected Growth):")
    print("The expected value grows linearly. The optimal liquidation strategy provides an expected return")
    print("that is *greater than* the expected value at the end of the period.")
    print("Lower Bound for Expected Value = Initial Investment * (1 + Growth Rate * Years)")
    print(f"E[Final Value] > {initial_investment:.1f} * (1 + {growth_rate} * {liquidation_year})")
    print(f"E[Final Value] > {expected_value_b_lower_bound:.4f}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Comparing Basket A's certain value ({value_a:.4f}) with the lower bound for Basket B's expected value ({expected_value_b_lower_bound:.4f}).")
    
    if value_a > expected_value_b_lower_bound:
        print("The value from compounding growth (Basket A) is significantly higher than the value from linear growth (Basket B).")
        print("Therefore, the strategy of choosing Basket A has a higher expected value.")
    else:
        # This case is not expected based on the math, but included for completeness.
        print("The expected value of Basket B is higher than the value of Basket A.")

analyze_baskets()