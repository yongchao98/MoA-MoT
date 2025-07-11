import sys
import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value of the minimum degree d_B,delta.

    Args:
        B (float): The upper bound of the interval [0, B].
        delta (float): The maximum allowed error.

    Returns:
        float: The asymptotic value A(B, delta).
    """
    if not B >= 1:
        raise ValueError("B must be >= 1")
    if not (0 < delta < 1):
        raise ValueError("delta must be in (0, 1)")

    # L is defined as log(delta^-1)
    L = math.log(1 / delta)

    # The asymptotic behavior is the sum of two terms.
    # Term 1 dominates when B is large compared to L.
    term1 = math.sqrt(B * L)
    
    # Term 2 dominates when L is large compared to B.
    # We use log(L + 1) to avoid math domain errors for L <= 1,
    # without changing the asymptotic behavior for large L.
    if L <= 0: # Should not happen given delta in (0,1)
        term2 = 0
    else:
        denominator_log = math.log(L + 1)
        if denominator_log == 0: # Handles L approaching 0
             # As L -> 0, L/log(L+1) -> 1. For simplicity, we can treat it as a small constant.
             # Or handle it based on the limit. lim_{x->0} x/log(x+1) = 1.
             term2 = 1.0
        else:
             term2 = L / denominator_log

    asymptotic_d = term1 + term2
    
    print(f"Given B = {B} and delta = {delta}")
    print(f"L = log(1/delta) = {L:.4f}")
    print("\nThe asymptotic formula for d_B,delta is: sqrt(B*L) + L/log(L)")
    print("We calculate it as: sqrt(B*L) + L/log(L+1) for stability.")
    print("\nBreaking down the calculation:")
    print(f"  Term 1 (sqrt(B*L)): sqrt({B} * {L:.4f}) = {term1:.4f}")
    print(f"  Term 2 (L/log(L+1)): {L:.4f} / log({L:.4f} + 1) = {term2:.4f}")
    print(f"  Result (Term 1 + Term 2): {term1:.4f} + {term2:.4f} = {asymptotic_d:.4f}")
    
    return asymptotic_d

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <B> <delta>")
        print("Example: python your_script_name.py 100 1e-10")
    else:
        try:
            B_val = float(sys.argv[1])
            delta_val = float(sys.argv[2])
            calculate_asymptotic_degree(B_val, delta_val)
        except ValueError as e:
            print(f"Error: {e}")
