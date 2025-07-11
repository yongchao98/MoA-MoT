import math

def analyze_baskets():
    """
    Calculates and compares the final value of Basket A and the expected 
    final value of Basket B over a 30-year period.
    """
    initial_value = 1.0
    rate = 0.10
    years = 30.0

    # --- Basket A Calculation ---
    # Basket A's value is determined by discrete annual compounding.
    value_a = initial_value * math.pow(1 + rate, years)

    # --- Basket B Calculation ---
    # Basket B's expected value is determined by continuous compounding of its growth trend.
    expected_value_b = initial_value * math.exp(rate * years)
    
    print("Comparing the value of Basket A with the expected value of Basket B after 30 years.")
    print("-" * 60)
    
    # Print the equation and result for Basket A
    print("Basket A is deterministic and grows with discrete annual compounding.")
    print(f"Final Value of Basket A = {initial_value} * (1 + {rate})^{int(years)}")
    print(f"Result: {value_a:.4f}")
    
    print("-" * 60)
    
    # Print the equation and result for Basket B
    print("Basket B is volatile and its expected value grows with continuous compounding.")
    print(f"Expected Final Value of Basket B = {initial_value} * e^({rate} * {int(years)})")
    print(f"Result: {expected_value_b:.4f}")
    
    print("-" * 60)
    
    # Compare the results and state the conclusion
    if expected_value_b > value_a:
        print("Conclusion: The expected value of Basket B is greater than the value of Basket A.")
    elif value_a > expected_value_b:
        print("Conclusion: The value of Basket A is greater than the expected value of Basket B.")
    else:
        print("Conclusion: The expected values of both baskets are equal.")

analyze_baskets()
