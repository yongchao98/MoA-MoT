import math

def analyze_generating_function(s, s_prime, kappa_s):
    """
    Analyzes the asymptotic behavior of the generating function H(s, s')
    for a planar billiard system.

    Args:
        s (float): The initial arc-length parameter.
        s_prime (float): The final arc-length parameter, close to s.
        kappa_s (float): The curvature of the boundary at s, κ(s).
    """
    if s == s_prime:
        print("H(s, s') = 0 when s = s'.")
        return

    delta_s = abs(s_prime - s)
    
    # The leading-order term is simply the arc-length separation.
    leading_term = delta_s
    
    # The curvature correction term.
    # Proportional to kappa^2 and delta_s^3.
    correction_term = (kappa_s**2 / 24) * (delta_s**3)
    
    # The asymptotic approximation of H(s, s')
    H_asymptotic = leading_term - correction_term
    
    print(f"Asymptotic Analysis of H(s, s') for s={s}, s'={s_prime}, κ(s)={kappa_s}")
    print("-" * 60)
    print("The formula for the generating function H(s,s') in the limit |s'-s| -> 0 is:")
    print("H(s,s') ≈ |s' - s| - (κ(s)² / 24) * |s' - s|³")
    print("\nCalculating each part:")
    
    # Print the full equation with numerical values substituted
    print(f"H({s}, {s_prime}) ≈ |{s_prime} - {s}| - ({kappa_s}² / 24) * |{s_prime} - {s}|³")
    
    print(f"           ≈ {delta_s:.6f} - ({kappa_s**2:.4f} / 24) * {delta_s**3:.6f}")
    
    # Print the calculated terms
    print(f"           ≈ {leading_term:.8f} - {correction_term:.8f}")
    
    # Print the final result
    print(f"           ≈ {H_asymptotic:.8f}")
    print("-" * 60)


if __name__ == '__main__':
    # Example usage:
    # Let's consider a point on the boundary s=1.0, and a nearby point s'=1.01
    # Assume the local curvature at s=1.0 is κ(1.0) = 0.5 (e.g., a circle of radius 2)
    s_val = 1.0
    s_prime_val = 1.01
    kappa_val = 0.5
    
    analyze_generating_function(s_val, s_prime_val, kappa_val)
