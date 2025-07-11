def print_potential_distribution_expression():
    """
    This function prints the derived mathematical expression for the EDL
    potential distribution psi(y) in a human-readable format.
    """
    
    # Define the components of the equation as strings for clear printing.
    # The expression is derived from the linearized Poisson-Boltzmann equation
    # with the boundary conditions psi(H/2) = 0 and psi(-H/2) = z_1*(1 + beta*k).
    
    numerator_zeta = "z_1 * (1 + beta * k)"
    denominator = "sinh(k * H)"
    sinh_term = "sinh(k * (H/2 - y))"
    
    # Construct the final expression string.
    # The numbers in the equation are 1 and 2.
    final_expression = f"psi(y) = (({numerator_zeta}) / ({denominator})) * {sinh_term}"
    
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    print(final_expression)
    print("\nWhere:")
    print("psi(y): Electrical potential at position y")
    print("z_1: Zeta potential parameter for the bottom surface (j=1)")
    print("beta: Slip length")
    print("k: Debye-Huckel parameter")
    print("H: Height of the microchannel")
    print("y: Vertical position in the channel (-H/2 <= y <= H/2)")
    print("sinh: Hyperbolic sine function")

# Execute the function to print the expression.
print_potential_distribution_expression()