import math

def analyze_generating_function(s_prime, s, kappa_s):
    """
    Calculates the asymptotic behavior of the billiard generating function H(s, s').

    This function characterizes the leading-order behavior of H(s, s')
    in the limit as the arc-length separation |s' - s| approaches zero.
    The analysis incorporates the local curvature kappa(s) of the boundary.

    The asymptotic formula is:
    H(s, s') ≈ |s' - s| - (κ(s)² / 24) * |s' - s|³

    Args:
        s_prime (float): The arc-length parameter of the destination point.
        s (float): The arc-length parameter of the initial point.
        kappa_s (float): The local curvature of the boundary at point s.
    """
    if s_prime == s:
        print("H(s, s') = 0 when s' = s.")
        return

    delta_s = s_prime - s
    abs_delta_s = abs(delta_s)

    # Term 1: The arc-length difference (leading term)
    term1_val = abs_delta_s

    # Term 2: The curvature correction
    # Coefficient part
    coeff_val = (kappa_s**2) / 24
    # Full term value
    term2_val = -coeff_val * (abs_delta_s**3)
    
    # Calculate the final result
    h_approx = term1_val + term2_val

    print(f"Asymptotic Analysis of H(s'={s_prime}, s={s}):")
    print("-" * 40)
    print("Formula: H(s, s') ≈ |s' - s| - (κ(s)² / 24) * |s' - s|³")
    print("\nCalculating each component:")

    # Print the equation with the specific numbers plugged in
    print(f"H ≈ |{s_prime} - {s}| - ({kappa_s}² / 24) * |{s_prime} - {s}|³")
    
    term1_str = f"{abs_delta_s:.5f}"
    kappa_sq_str = f"{kappa_s**2:.4f}"
    delta_s_cubed_str = f"{abs_delta_s**3:.7f}"
    
    print(f"H ≈ {term1_str} - ({kappa_sq_str} / 24) * {delta_s_cubed_str}")

    term2_calc_str = f"(-{coeff_val:.6f}) * {delta_s_cubed_str}"
    print(f"H ≈ {term1_str} + ({term2_calc_str})")

    # Print the final numerical result
    print("-" * 40)
    print(f"Result: H({s_prime}, {s}) ≈ {h_approx:.8f}")

# --- Example Usage ---
# Consider a point s=1.0 on a boundary with local curvature kappa(s)=2.0.
# We want to find the generating function value for a nearby point s'=1.05.
s_val = 1.0
s_prime_val = 1.05
kappa_at_s = 2.0  # Example curvature, e.g., for a circle of radius 0.5

analyze_generating_function(s_prime_val, s_val, kappa_at_s)
