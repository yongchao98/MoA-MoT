import math

def analyze_generating_function():
    """
    This script presents the result of the asymptotic analysis of the
    billiard generating function H(s, s') in the limit of small
    separation |s' - s|.
    """

    # The analysis provides the coefficients and powers for the expansion.
    # The first term represents the distance along the tangent (a straight line).
    coeff_term1 = 1
    power_term1 = 1

    # The second term is the first-order correction due to the boundary's
    # curvature, kappa(s). The derivation shows its coefficient is 1/24.
    coeff_term2_num = 1
    coeff_term2_den = 24

    # The correction term's dependence on the separation |s' - s| is cubic.
    power_term2 = 3

    # The order of the next term in the expansion.
    # This is often referred to as the Big-O term.
    power_big_o = 4

    # Construct and print the final equation for H(s, s').
    # The equation shows how the simple distance is corrected by the
    # geometric curvature of the billiard boundary.
    print("The leading-order asymptotic expansion of the generating function H(s, s') is:")

    equation = (
        f"H(s, s')  ≈  {coeff_term1} * |s' - s|^{power_term1}  -  "
        f"({coeff_term2_num}/{coeff_term2_den}) * κ(s)² * |s' - s|^{power_term2}  +  "
        f"O(|s' - s|^{power_big_o})"
    )

    print(equation)
    print("\nWhere:")
    print("  s, s'   : arc-length parameters of two consecutive collision points.")
    print("  κ(s)    : the local curvature of the boundary at point s.")
    print("  O(...)  : represents the higher-order terms in the expansion.")

if __name__ == "__main__":
    analyze_generating_function()
