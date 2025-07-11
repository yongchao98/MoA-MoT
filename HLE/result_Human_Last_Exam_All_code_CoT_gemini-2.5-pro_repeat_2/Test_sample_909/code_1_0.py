def solve_and_print_electric_field():
    """
    This function prints the derived formulas for the electric field in each region.
    The derivation is based on solving Laplace's equation with the given boundary conditions.
    """

    # The problem is symbolic, so we will construct the string representation of the formulas.
    # The variables are:
    # V0: Applied DC voltage
    # sigma1, sigma2: Ohmic conductivities of the two regions
    # r: radial distance from the center
    # i_phi: unit vector in the azimuthal direction
    # pi: mathematical constant

    # Expression for the electric field in Region 1 (0 < phi < pi/2)
    E1_expression = "E_1 = (2 * sigma_2 * V0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"

    # Expression for the electric field in Region 2 (pi/2 < phi < pi)
    E2_expression = "E_2 = (2 * sigma_1 * V0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"

    print("The electric field in Region 1 is:")
    print(E1_expression)
    print("\nThe electric field in Region 2 is:")
    print(E2_expression)

# Execute the function to print the results
solve_and_print_electric_field()