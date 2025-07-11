import sympy as sp

def solve_edl_potential():
    """
    This function defines and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y)
    in a parallel-plate microchannel under the given conditions.
    """
    # Define the symbols used in the equation
    # y: spatial coordinate across the channel height
    # H: height of the microchannel
    # k: Debye-Huckel parameter
    # beta: slip length
    # z1: zeta potential parameter for the bottom surface
    y, H, k, beta, z1 = sp.symbols('y H k beta z_1')
    
    # psi(y) represents the EDL potential as a function of y
    psi = sp.Function('psi')(y)

    # The slip-dependent zeta potential at the bottom surface (y = -H/2) is z_a1.
    # The potential at the top surface (y = H/2) is 0 since z_2 = 0.
    z_a1 = z1 * (1 + beta * k)

    # The final expression for the EDL potential distribution is derived by solving
    # the linearized Poisson-Boltzmann equation d^2(psi)/dy^2 = k^2 * psi
    # with the boundary conditions: psi(H/2) = 0 and psi(-H/2) = z_a1.
    # The resulting expression is:
    # psi(y) = z_a1 * sinh(k*(H/2 - y)) / sinh(k*H)
    
    edl_potential_expression = z_a1 * sp.sinh(k * (H / 2 - y)) / sp.sinh(k * H)

    # Create a symbolic equation psi(y) = expression for clear display
    final_equation = sp.Eq(psi, edl_potential_expression)

    # Print the final expression for the user
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # sp.pretty_print provides a nicely formatted output in the console
    sp.pretty_print(final_equation)

if __name__ == '__main__':
    solve_edl_potential()