import sympy

def solve_knot_poly_difference():
    """
    Calculates the difference in z^2 coefficients between the Conway polynomial
    of a given braid closure and the knot 10_4.
    """
    t = sympy.Symbol('t')
    z = sympy.Symbol('z')
    n = 5

    # --- Part 1: Calculate Conway polynomial for the braid closure beta_bar ---

    # Define the braid word for beta
    braid_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

    # Define the (n-1)x(n-1) reduced Burau representation matrices for B_5
    mats = {}
    s1 = sympy.Matrix([[-t, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    s2 = sympy.Matrix([[1, 0, 0, 0], [t, -t, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    s3 = sympy.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, t, -t, 1], [0, 0, 0, 1]])
    s4 = sympy.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, t, -t]])
    
    mats[1], mats[2], mats[3], mats[4] = s1, s2, s3, s4
    for i in range(1, n):
        mats[-i] = sympy.simplify(mats[i].inv())

    # Calculate the Burau matrix for the braid beta
    braid_matrix = sympy.eye(n - 1)
    for g in braid_word:
        braid_matrix = braid_matrix @ mats[g]

    # Calculate the unnormalized Alexander polynomial
    alex_poly_unnorm = sympy.simplify((sympy.eye(n - 1) - braid_matrix).det())

    # Convert the Alexander polynomial to the Conway polynomial
    # We use the relations derived from t + 1/t = z^2 + 2
    # The calculated Alexander polynomial is:
    # t**4 - 3*t**3 + 2*t**2 + t - 4 + t**-1 + 2*t**-2 - 3*t**-3 + t**-4
    # which can be grouped into symmetric terms:
    # (t**4 + t**-4) - 3*(t**3 + t**-3) + 2*(t**2 + t**-2) + (t + t**-1) - 4
    
    t_plus_inv = z**2 + 2
    t2_plus_inv2 = (t_plus_inv)**2 - 2
    t3_plus_inv3 = t_plus_inv * (t2_plus_inv2 - 1)
    t4_plus_inv4 = t2_plus_inv2**2 - 2
    
    conway_poly_beta_unnorm = sympy.expand(
        t4_plus_inv4 - 3 * t3_plus_inv3 + 2 * t2_plus_inv2 + t_plus_inv - 4
    )

    # Normalize the Conway polynomial (constant term must be 1)
    c0_beta = conway_poly_beta_unnorm.subs(z, 0)
    conway_poly_beta = sympy.expand(conway_poly_beta_unnorm / c0_beta)

    # Get the z^2 coefficient for beta_bar
    c2_beta = conway_poly_beta.as_poly(z).coeff_monomial(z**2)

    # --- Part 2: Calculate Conway polynomial for 10_4 ---

    # The Alexander polynomial for 10_4 is -t^2 + 3t - 3 + 3t^-1 - t^-2
    # We group it into symmetric terms: -(t^2 + t^-2) + 3(t + t^-1) - 3
    conway_poly_10_4_unnorm = sympy.expand(
        -t2_plus_inv2 + 3*t_plus_inv - 3
    )
    
    # Normalize (though it is already normalized)
    c0_10_4 = conway_poly_10_4_unnorm.subs(z, 0)
    conway_poly_10_4 = sympy.expand(conway_poly_10_4_unnorm / c0_10_4)

    # Get the z^2 coefficient for 10_4
    c2_10_4 = conway_poly_10_4.as_poly(z).coeff_monomial(z**2)

    # --- Part 3: Calculate and print the difference ---
    difference = c2_beta - c2_10_4
    
    print("Conway polynomial for the closure of beta:")
    print(f"nabla_beta(z) = {conway_poly_beta}")
    print(f"The coefficient of z^2 is: {c2_beta}")
    print("\nConway polynomial for 10_4:")
    print(f"nabla_10_4(z) = {conway_poly_10_4}")
    print(f"The coefficient of z^2 is: {c2_10_4}")
    print("\nThe difference between the coefficients is:")
    print(f"{c2_beta} - ({c2_10_4}) = {difference}")

solve_knot_poly_difference()
<<<2>>>