def solve_magnetization_curve():
    """
    This function provides the analytical expression for the initial magnetization
    curve of a superconducting bar as described in the problem.
    It prints the derivation steps and the final formula.
    """

    print("Derivation for the initial magnetization curve M(H):")
    print("-----------------------------------------------------")
    print("The problem considers a superconducting bar with a very wide, thin rectangular")
    print("cross-section (-a <= x <= a, -b <= y <= b with b << a) in a parallel")
    print("applied magnetic field H along the y-axis.")
    print("\nStep 1: The system is modeled as an infinite slab of width 2a in a parallel field H.")
    print("\nStep 2: The Bean critical-state model is used, with a constant critical current density Jc.")
    print("         This model describes the state for H increasing from a zero-field-cooled state.")
    print("\nStep 3: The magnetization M is calculated using the standard definition for bulk samples:")
    print("         M = <H_internal> - H_applied")
    print("         This is valid for applied fields H below the full penetration field, H_p = a * Jc.")

    print("\nFinal Analytical Expression:")
    print("-----------------------------")
    print("The resulting equation for the magnetization M as a function of the applied field H is:")
    # M(H) = H^2 / (2*a*Jc) - H
    
    # We will print the equation piece by piece to highlight the numbers.
    # The coefficients in the polynomial form are 1/(2*a*Jc) and -1.
    
    numerator_coeff = 1
    denominator_coeff = 2
    linear_coeff = -1

    print("\nThe equation is parabolic:")
    print(f"M(H) = (H^2 / ({denominator_coeff} * a * J_c)) - H")
    
    print("\nLet's break down the components and numerical coefficients:")
    print("Term 1 (Quadratic term): H^2 / (2 * a * J_c)")
    print(f"  - The numerical coefficient in the denominator is: {denominator_coeff}")
    
    print("Term 2 (Linear term): -H")
    print(f"  - The numerical coefficient is: {linear_coeff}")
    
    print("\nWhere:")
    print("  M(H): Magnetization")
    print("  H: Applied magnetic field")
    print("  a: Half-width of the superconducting bar")
    print("  J_c: Critical current density")


solve_magnetization_curve()