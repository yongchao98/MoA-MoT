import sympy

def display_generating_function_asymptotics():
    """
    This function derives and displays the asymptotic expansion of the
    billiard generating function H(s,s') for s' -> s.
    The final expression is printed to the console.
    """

    # Define symbolic representations for the mathematical terms.
    # These are for display purposes.
    H_ss_prime = "H(s, s')"
    s_diff_abs = "|s' - s|"
    kappa_s_sq = "κ(s)²"

    # The coefficients and powers from the derivation.
    leading_term_coeff = 1
    correction_term_coeff_num = 1
    correction_term_coeff_den = 24
    s_diff_power_in_correction = 3
    error_order_power = 5

    # Build the equation string piece by piece to satisfy the prompt's
    # requirement to "output each number in the final equation".
    
    # Term 1: The trivial leading term, which is the arc length.
    # Note: its coefficient is 1.
    term1 = f"{s_diff_abs}"

    # Term 2: The first non-trivial correction due to curvature.
    # Numbers involved: 1, 24, 3.
    term2 = f"({correction_term_coeff_num}/{correction_term_coeff_den}) * {kappa_s_sq} * {s_diff_abs}**{s_diff_power_in_correction}"

    # Term 3: The order of the error term.
    # Number involved: 5.
    error_term = f"O({s_diff_abs}**{error_order_power})"
    
    # Construct the final equation.
    equation = f"{H_ss_prime} ≈ {term1} - {term2} + {error_term}"
    
    print("Within the framework of planar Birkhoff billiard dynamics, the leading-order asymptotic expansion of the generating function H(s, s') in the limit |s'-s| -> 0 is given by:")
    print("-" * 80)
    print(equation)
    print("-" * 80)
    print("\nThis formula shows that the primary deviation from simple arc-length distance is a third-order term that depends on the square of the local boundary curvature κ(s).")
    
    print("\nBreakdown of numerical components in the equation:")
    print(f"- Leading Term Coefficient: {leading_term_coeff}")
    print(f"- Correction Term Coefficient: -{correction_term_coeff_num}/{correction_term_coeff_den}")
    print(f"- Power in Correction Term: {s_diff_power_in_correction}")
    print(f"- Order of Higher Terms: {error_order_power}")


if __name__ == "__main__":
    display_generating_function_asymptotics()
