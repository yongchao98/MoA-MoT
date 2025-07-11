import sympy
from sympy import symbols, Function, sinh, Eq, pprint

def solve_edl_potential():
    """
    This function derives and prints the expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) in a parallel-plate microchannel with asymmetric
    zeta potentials, where one plate has zero potential.
    """
    # Define the symbolic variables and the potential function psi(y)
    y, z_1, beta, k, H = symbols('y zeta_1 beta k H', real=True, positive=True)
    psi = Function('psi')(y)

    # From the problem statement, the slip-dependent zeta potential at the bottom surface (j=1) is:
    # z_aj = z_j(1 + beta*k)
    # The potential at the top surface (j=2) is zero (z_2 = 0).
    z_a1 = z_1 * (1 + beta * k)

    # The derived solution for psi(y) with the coordinate system centered in the channel
    # (y ranging from -H/2 to H/2) is:
    # psi(y) = z_a1 * sinh(k*(H/2 - y)) / sinh(k*H)
    # This form is obtained by solving the linearized Poisson-Boltzmann equation d^2(psi)/dy^2 = k^2*psi
    # with boundary conditions: psi(-H/2) = z_a1 and psi(H/2) = 0.
    
    # Construct the right-hand side of the equation
    numerator = z_a1 * sinh(k * (H/2 - y))
    denominator = sinh(k * H)
    rhs = numerator / denominator

    # Create the full equation object for pretty printing
    final_equation = Eq(psi, rhs)

    # Print the final expression in a human-readable format.
    # The numbers 1 and 2 are present in the terms (1 + beta*k) and H/2.
    print("The final expression for the Electrical Double-Layer (EDL) potential distribution psi(y) is:")
    pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    solve_edl_potential()