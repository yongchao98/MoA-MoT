def print_magnetic_field_solution():
    """
    This function prints the derived magnetic field H(r, theta) for both regions.
    The formula is presented in a clear, readable format corresponding to the final answer.
    """
    
    # Header for the solution
    print("The magnetic field H(r, θ) is given by:")
    print("-" * 40)
    
    # Case 1: Inside the sphere (0 < r < R)
    # H_in = (2 * mu_0 / mu) * (K_0 / (1 + (2 * mu_0 / mu))) * z_hat
    print("For the region inside the sphere (0 < r < R):")
    
    # Constructing the string for the inside field
    # Numerator part: 2 * mu_0 * K_0
    numerator_in = "2 \u03bc\u2080 K\u2080"
    # Denominator part: mu * (1 + 2 * mu_0 / mu)
    denominator_in = "\u03bc (1 + 2 \u03bc\u2080/\u03bc)"
    # z_hat vector
    z_hat = "z\u0302"
    
    print(f"  \u20d7")
    print(f"  H_in  =  ({numerator_in} / {denominator_in}) {z_hat}")
    print("\nWhich can be rewritten as:\n")
    print(f"  \u20d7")
    print(f"  H_in  =  (2 \u03bc\u2080/\u03bc) * (K\u2080 / (1 + 2 \u03bc\u2080/\u03bc)) {z_hat}")
    print("\n")
    
    # Case 2: Outside the sphere (r > R)
    # H_out = (K_0 / (1 + (2 * mu_0 / mu))) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    print("For the region outside the sphere (r > R):")

    # Constructing the string for the outside field
    # First term: K_0 / (1 + 2 * mu_0 / mu)
    term1_out = "K\u2080 / (1 + 2 \u03bc\u2080/\u03bc)"
    # Second term: R^3 / r^3
    term2_out = "R\u00b3 / r\u00b3"
    # Vector part: (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    vector_part_out = "(2 cos(\u03b8) r\u0302 + sin(\u03b8) \u03b8\u0302)"

    print(f"  \u20d7")
    print(f"  H_out =  ({term1_out}) * ({term2_out}) * {vector_part_out}")
    print("-" * 40)

# Execute the function to print the solution
print_magnetic_field_solution()