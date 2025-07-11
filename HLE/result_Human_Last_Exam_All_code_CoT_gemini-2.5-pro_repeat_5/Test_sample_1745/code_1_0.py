import sympy

def solve_edl_potential():
    """
    This function derives and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y)
    in a parallel-plate microchannel under the specified conditions.
    """
    # Define the symbolic variables
    # psi: The electrical potential as a function of y
    # y: The spatial coordinate across the channel height
    # z_1: The zeta potential parameter for the bottom surface
    # beta: The slip length
    # k: The Debye-Huckel parameter
    # H: The height of the microchannel
    psi = sympy.Function('psi')
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H')

    # The problem asks for the expression for the EDL potential distribution.
    # Based on the Poisson-Boltzmann equation and the boundary conditions:
    # 1. psi(0) = z_1 * (1 + beta * k)
    # 2. psi(H) = 0
    # The solution is derived as:
    # psi(y) = z_1*(1 + beta*k) * sinh(k*(H-y)) / sinh(k*H)

    # Construct the final expression
    numerator = sympy.sinh(k * (H - y))
    denominator = sympy.sinh(k * H)
    prefactor = z_1 * (1 + beta * k)

    rhs_expression = prefactor * (numerator / denominator)

    # Create the equation psi(y) = expression
    final_equation = sympy.Eq(psi(y), rhs_expression)

    # Print the final equation in a readable format
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    sympy.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    solve_edl_potential()
