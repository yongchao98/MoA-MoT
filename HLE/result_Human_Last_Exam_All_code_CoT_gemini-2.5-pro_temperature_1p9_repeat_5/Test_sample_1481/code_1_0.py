def print_generating_function_asymptotic():
    """
    This function prints the leading-order asymptotic expansion of the
    billiard generating function H(s,s') for |s' - s| -> 0.

    The expansion reveals the influence of the boundary's local
    curvature, kappa(s), on the system's dynamics.
    """
    
    # Define symbolic representations for the terms in the equation.
    gen_func = "H(s, s')"
    separation = "|s' - s|"
    curvature_term = "κ(s)²"
    
    # Define the coefficients from the derivation.
    # Note: Python's float representation would be 0.041666...
    # The fractional form is more precise and illustrative.
    coeff_numerator = 1
    coeff_denominator = 24
    
    # Construct the final equation string using an f-string.
    # The equation shows the generating function is approximated by the arc-length
    # separation minus a term dependent on the cube of the separation and
    # the square of the curvature.
    asymptotic_expansion = (
        f"{gen_func} ≈ {separation} - "
        f"({coeff_numerator}/{coeff_denominator}) * {curvature_term} * {separation}³ + "
        f"O({separation}⁵)"
    )

    print("The asymptotic expansion for the billiard generating function H(s,s') is:")
    print(asymptotic_expansion)

if __name__ == "__main__":
    print_generating_function_asymptotic()
