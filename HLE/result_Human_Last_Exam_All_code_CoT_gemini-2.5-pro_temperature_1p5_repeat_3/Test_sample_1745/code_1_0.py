import sympy

def solve_potential_distribution():
    """
    This function defines the symbolic expression for the electrical double-layer
    potential distribution psi(y) and prints it.
    """
    # Define the mathematical symbols involved in the equation.
    # y: the spatial coordinate in the direction perpendicular to the channel walls
    # z_1: the zeta potential of the bottom surface
    # beta: the slip length
    # k: the Debye-Huckel parameter
    # H: the height of the microchannel
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H')

    # The expression for the slip-dependent zeta potential at the bottom surface (j=1)
    # is z_a1 = z_1 * (1 + beta * k). The number '1' is included here.
    slip_dependent_zeta = z_1 * (1 + beta * k)

    # The final expression for the EDL potential distribution psi(y) is derived from
    # the linearized Poisson-Boltzmann equation with the given boundary conditions.
    # psi(y) = zeta_at_wall * sinh(k*(H-y)) / sinh(k*H)
    potential_distribution = slip_dependent_zeta * (sympy.sinh(k * (H - y)) / sympy.sinh(k * H))

    # Create a symbolic function for psi(y) to create a formal equation.
    psi_of_y = sympy.Function('psi')(y)

    # Form the final equation object for printing.
    final_equation = sympy.Eq(psi_of_y, potential_distribution)

    # Print the final equation in a readable format.
    print("The final expression for the Electrical double-layer potential distribution is:")
    # Using sympy.pprint to render the mathematical formula.
    # It ensures all symbols and numbers like '1' are correctly displayed.
    sympy.pprint(final_equation, use_unicode=True)

if __name__ == "__main__":
    solve_potential_distribution()