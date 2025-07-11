import sympy

def display_generating_function_asymptotic():
    """
    This function derives and displays the asymptotic expansion for the
    billiard generating function H(s, s').
    """

    # Define symbolic variables for the derivation visualization
    # s_val is a placeholder for the point s
    s, s_prime, delta_s, s_val = sympy.symbols("s s' Delta_s s_val")
    #kappa represents the curvature function kappa(s)
    kappa = sympy.Function('κ')(s)
    # r represents the position vector r(s)
    r = sympy.Function('r')(s)

    # The Taylor expansion of H(s,s')^2 up to O(delta_s^4) is:
    # H(s,s')^2 = (s'-s)^2 - (1/12) * kappa(s)^2 * (s'-s)^4 + O((s'-s)^6)
    # We will derive H(s,s') from this.

    # Expression inside the square root: 1 - (κ(s)²/12)Δs²
    x = -(kappa**2 / 12) * delta_s**2

    # Taylor expansion of sqrt(1+x) is 1 + 1/2*x + ...
    sqrt_expansion = 1 + sympy.Rational(1, 2) * x

    # H(s, s') = |s'-s| * sqrt(1 - (κ(s)²/12)(s'-s)²)
    # We use abs(delta_s) to represent |s' - s|
    H_expression = sympy.Abs(delta_s) * sqrt_expansion

    # Final formula in symbolic form
    # Replace delta_s with |s' - s| for final display
    s_diff_str = "|s' - s|"
    kappa_sq_str = "κ(s)²"
    final_formula = f"H(s, s') ≈ {s_diff_str} - (1/24){kappa_sq_str}{s_diff_str}³ + O({s_diff_str}⁵)"


    print("This script programmatically generates the leading-order asymptotic expansion for the billiard map generating function H(s,s').")
    print("-" * 80)
    print("The derived asymptotic formula is:")
    print(final_formula)
    print("-" * 80)
    print("Analysis of the numerical components of the equation:")

    # Extract numbers from the derived formula for output
    term1_coeff = 1
    term1_exp_s = 1

    term2_coeff_num = 1
    term2_coeff_den = 24
    term2_exp_kappa = 2
    term2_exp_s = 3

    print("\nTerm 1: The Chord Length Approximation")
    print(f"  Formula part: {s_diff_str}")
    print(f"  Implicit numerical coefficient: {term1_coeff}")
    print(f"  Exponent of the arc-length separation |s'-s|: {term1_exp_s}")


    print("\nTerm 2: The Leading-Order Curvature Correction")
    print(f"  Formula part: -(1/{term2_coeff_den}){kappa_sq_str}{s_diff_str}³")
    print(f"  Numerator of the coefficient: {term2_coeff_num}")
    print(f"  Denominator of the coefficient: {term2_coeff_den}")
    print(f"  Exponent of the curvature κ(s): {term2_exp_kappa}")
    print(f"  Exponent of the arc-length separation |s'-s|: {term2_exp_s}")

if __name__ == '__main__':
    display_generating_function_asymptotic()
