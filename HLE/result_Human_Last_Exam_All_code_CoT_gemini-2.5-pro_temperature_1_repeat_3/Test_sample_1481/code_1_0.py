def display_generating_function_asymptotic():
    """
    This script presents the result of the asymptotic analysis of the
    generating function H(s, s') for a planar Birkhoff billiard system
    in the limit as the arc-length separation |s' - s| approaches zero.
    """

    # Define the symbolic components of the equation for clarity
    H_ss_prime = "H(s, s')"
    s_prime_minus_s = "|s' - s|"
    kappa_s_sq = "κ(s)²"
    s_prime_minus_s_cubed = "|s' - s|³"
    big_O_term = "O(|s' - s|⁵)"

    # Coefficients derived from the analysis
    coeff_first_term = 1
    coeff_second_term_numerator = 1
    coeff_second_term_denominator = 24

    # Print the explanation and the final formula
    print("This analysis characterizes the generating function H(s,s'), which equals the")
    print("chord length between two points q(s) and q(s') on the billiard boundary.")
    print("In the asymptotic limit where the separation |s' - s| is very small,")
    print("the function can be expanded in terms of the local boundary curvature κ(s).")
    print("\nThe leading-order behavior is given by the following equation:")
    print("-" * 70)

    # H(s, s') = 1 * |s' - s| - (1/24) * κ(s)² * |s' - s|³ + O(|s' - s|⁵)
    print(f"{H_ss_prime} = {coeff_first_term} * {s_prime_minus_s} - "
          f"({coeff_second_term_numerator}/{coeff_second_term_denominator}) * "
          f"{kappa_s_sq} * {s_prime_minus_s_cubed} + {big_O_term}")

    print("-" * 70)

if __name__ == "__main__":
    display_generating_function_asymptotic()