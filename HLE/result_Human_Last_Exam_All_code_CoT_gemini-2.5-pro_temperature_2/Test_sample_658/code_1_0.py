import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value of the minimum degree d_B_delta
    for polynomial approximation of exp(-x) on [0, B] with error delta.

    Args:
        B (float): The upper bound of the interval [0, B].
        delta (float): The maximum allowed error.
    """
    if B <= 0 or delta <= 0 or delta >= 1:
        print("Please provide B > 0 and 0 < delta < 1.")
        return

    # L is defined as log(1/delta)
    L = math.log(1/delta)

    # First term in the asymptotic formula
    term1 = math.sqrt(B * L)

    # Second term in the asymptotic formula
    # We use log(L/B + 2) to ensure the argument of the log is > 1
    # as L/B can be small. For asymptotic analysis, the constant is not critical.
    log_arg = L/B + 2
    term2 = L / math.log(log_arg)
    
    # The asymptotic degree is the sum of the two terms
    d_asymptotic = term1 + term2

    print(f"Given B = {B} and delta = {delta}")
    print(f"Calculated L = log(1/delta) = {L:.4f}")
    print("The asymptotic value is composed of two terms:")
    print(f"Term 1 (sqrt(B*L)): {term1:.4f}")
    print(f"Term 2 (L/log(L/B+2)): {term2:.4f}")
    print(f"Final asymptotic value for d_{B,delta}: {d_asymptotic:.4f}")

# Example Usage:
# Regime 1: B >> L
print("--- Example 1: B is large compared to L ---")
calculate_asymptotic_degree(B=100, delta=1e-10)

print("\n" + "="*40 + "\n")

# Regime 2: L >> B
print("--- Example 2: L is large compared to B ---")
calculate_asymptotic_degree(B=10, delta=1e-100)