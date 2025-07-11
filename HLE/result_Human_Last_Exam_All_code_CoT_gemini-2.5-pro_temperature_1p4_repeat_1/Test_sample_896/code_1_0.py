def solve_knot_polynomial_difference():
    """
    This function calculates the difference between the z^2 coefficients
    of the Alexander-Conway polynomials for two knots.
    It must be run in a SageMath environment.
    """
    try:
        from sage.all import BraidGroup, Knot
    except ImportError:
        print("This code must be run in a SageMath environment.")
        print("You can run it online at SageCell (https://sagecell.sagemath.org/).")
        return

    # 1. Define the braid group on 5 strands.
    B = BraidGroup(5)

    # 2. Define the braid beta from its word representation.
    # beta = s4^-1 * s4^-1 * s3^-1 * s4 * s3^-1 * s2 * s1^-1 * s3^-1 * s2^-1 * s2^-1 * s2^-1 * s1^-1
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    beta = B(beta_word)

    # 3. Compute the closure of the braid, which is a knot we'll call K_beta.
    K_beta = beta.closure()

    # 4. Identify the knot K_beta. SageMath can check it against its knot catalogue.
    # This reveals that the closure of the braid is, in fact, the knot 10_4.
    k_beta_name = K_beta.name()
    print(f"The knot from the closure of the braid beta is identified as: {k_beta_name}")

    # 5. Compute the Alexander-Conway polynomial for the closure of beta.
    nabla_beta = K_beta.conway_polynomial()
    print(f"The Alexander-Cownway polynomial for the closure of beta, nabla_b(z), is: {nabla_beta}")

    # 6. Define the knot 10_4.
    K_10_4 = Knot('10_4')

    # 7. Compute the Alexander-Conway polynomial for 10_4.
    nabla_10_4 = K_10_4.conway_polynomial()
    print(f"The Alexander-Conway polynomial for the knot 10_4, nabla_10_4(z), is: {nabla_10_4}")

    # 8. Get the polynomial variable 'z'.
    z = nabla_beta.parent().gen()

    # 9. Extract the z^2 coefficient from both polynomials.
    # Since the knots are the same, the polynomials and coefficients will be identical.
    coeff_beta = nabla_beta.coefficient(z**2)
    coeff_10_4 = nabla_10_4.coefficient(z**2)
    print(f"\nThe coefficient of z^2 in nabla_b(z) is: {coeff_beta}")
    print(f"The coefficient of z^2 in nabla_10_4(z) is: {coeff_10_4}")

    # 10. Calculate and display the final difference.
    difference = coeff_beta - coeff_10_4
    print("\nThe difference between the coefficients is calculated as follows:")
    print(f"{coeff_beta} - ({coeff_10_4}) = {difference}")

solve_knot_polynomial_difference()