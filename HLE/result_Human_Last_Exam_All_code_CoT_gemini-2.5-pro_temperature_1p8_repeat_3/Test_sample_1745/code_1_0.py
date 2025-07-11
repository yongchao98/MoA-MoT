def solve_edl_potential_distribution():
    """
    This function prints the derived analytical expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y)
    for the given boundary conditions.
    """

    # Define the components of the equation as strings
    lhs = "psi(y)"
    z_a1 = "z_1 * (1 + beta * k)"
    numerator = "sinh(k * (H/2 - y))"
    denominator = "sinh(k * H)"

    # Construct the final equation string
    # The equation is: psi(y) = z_1 * (1 + beta * k) * sinh(k * (H/2 - y)) / sinh(k * H)
    # The numbers in this equation are 1 (from z_1), 1 (from 1 + beta*k), and 2 (from H/2)
    final_expression = f"{lhs} = {z_a1} * ( {numerator} / {denominator} )"

    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    print(final_expression)

solve_edl_potential_distribution()