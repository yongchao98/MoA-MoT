def display_billiard_generating_function_asymptotics():
    """
    This function programmatically constructs and prints the final derived
    asymptotic formula for the planar billiard generating function H(s, s').
    The code explicitly outputs each number involved in the final equation.
    """

    # These are the numerical constants and powers derived from the
    # asymptotic analysis of the generating function.

    # Leading term: C1 * |s'-s|^p1
    coeff_leading = 1
    power_leading = 1

    # Curvature correction term: C2 * κ(s)^p_kappa * |s'-s|^p_ds
    coeff_correction_num = 1
    coeff_correction_den = 24
    power_curvature = 2
    power_separation = 3

    # Higher order term power
    power_big_o = 5
    
    # Using ds as a shorthand for the arc-length separation |s'-s|
    # and k(s) as a shorthand for curvature κ(s).
    
    # We construct the final equation string, outputting each number
    # directly into the formula.
    # Note: It's conventional to omit |s'-s|^1 and just write |s'-s|.
    
    final_equation = (
        f"H(s, s') ≈ {coeff_leading}*|s' - s| "
        f"- ({coeff_correction_num}/{coeff_correction_den})*κ(s)^{power_curvature}*|s' - s|^{power_separation} "
        f"+ O(|s' - s|^{power_big_o})"
    )

    print("The asymptotic analysis of the planar billiard generating function H(s, s')")
    print("yields the following leading-order behavior as |s' - s| → 0:")
    print("-" * 70)
    print(final_equation)
    print("-" * 70)

display_billiard_generating_function_asymptotics()