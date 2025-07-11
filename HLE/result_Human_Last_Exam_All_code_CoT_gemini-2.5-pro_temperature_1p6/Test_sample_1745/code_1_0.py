import sympy as sp

def solve_potential_distribution():
    """
    This function derives and prints the symbolic expression for the EDL potential
    distribution in a parallel-plate microchannel based on the given problem.
    """
    # Define the symbols used in the equation.
    # y: vertical coordinate
    # H: height of the microchannel
    # k: Debye-Huckel parameter
    # beta: slip length
    # zeta_1: zeta potential of the bottom surface
    # psi(y): the EDL potential as a function of y
    y, H, k, beta = sp.symbols('y H k beta')
    zeta_1 = sp.Symbol('zeta_1')
    psi = sp.Function('psi')

    # Define the slip-dependant zeta potential at the bottom wall (j=1).
    # The potential at the top wall (j=2) is zero since zeta_2 = 0.
    zeta_a1 = zeta_1 * (1 + beta * k)

    # The solution to the linearized Poisson-Boltzmann equation d²ψ/dy² = k²ψ
    # with boundary conditions ψ(0) = ζ_a1 and ψ(H) = 0 is:
    # ψ(y) = ζ_a1 * sinh(k(H-y)) / sinh(kH)
    rhs = zeta_a1 * (sp.sinh(k * (H - y)) / sp.sinh(k * H))

    # Create the final equation for the EDL potential distribution.
    # The equation includes all the specified parameters: ζ₁, β, k, H, and the variable y.
    final_equation = sp.Eq(psi(y), rhs)

    # Print the final equation in a nicely formatted way.
    print("The expression for the Electrical double-layer potential distribution ψ(y) is:")
    sp.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    solve_potential_distribution()