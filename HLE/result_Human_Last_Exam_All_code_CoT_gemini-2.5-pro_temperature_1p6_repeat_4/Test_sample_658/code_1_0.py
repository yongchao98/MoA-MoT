import math

def calculate_asymptotic_d(B, delta):
    """
    Calculates the asymptotic value of d_B_delta.

    Args:
        B (float): The upper bound of the interval [0, B]. Must be >= 1.
        delta (float): The desired error bound. Must be in (0, 1).

    Returns:
        float: The asymptotic value of d_B_delta.
    """
    if B < 1 or not (0 < delta < 1):
        print("Invalid input: B must be >= 1 and delta must be in (0, 1).")
        return None
    
    # L is defined as log(1/delta)
    L = math.log(1 / delta)

    # To avoid issues with log(0) or division by zero when B is very large or L/B is small,
    # we add a constant, e, to the argument of the logarithm.
    # The choice of 'e' ensures log(L/B + e) is at least 1 when L/B approaches 0.
    log_regularizer = math.e
    
    # First term of the asymptotic formula, dominant when B >> L
    term1 = math.sqrt(B * L)
    
    # Second term, dominant when L >> B
    # It's L / log(L/B + e)
    term2 = L / math.log(L / B + log_regularizer)
    
    # The asymptotic value is the sum of these two terms.
    result = term1 + term2
    
    print(f"Given B = {B} and delta = {delta}")
    print(f"L = log(1/delta) = {L:.4f}")
    print("\nThe asymptotic formula for d_B_delta is composed of two terms:")
    print(f"Term 1 (sqrt(B*L)): {term1:.4f}")
    print(f"Term 2 (L/log(L/B + e)): {term2:.4f}")
    print(f"\nAsymptotic d_B_delta (Term 1 + Term 2): {result:.4f}")

    return result

# Example 1: B >> L regime
print("--- Example 1: B is large compared to L ---")
B1 = 1000.0
delta1 = 1e-4 # This gives L = log(10000) approx 9.21
calculate_asymptotic_d(B1, delta1)

print("\n" + "="*50 + "\n")

# Example 2: L >> B regime
print("--- Example 2: L is large compared to B ---")
B2 = 2.0
delta2 = 1e-50 # This gives L = log(1e50) approx 115.13
calculate_asymptotic_d(B2, delta2)
