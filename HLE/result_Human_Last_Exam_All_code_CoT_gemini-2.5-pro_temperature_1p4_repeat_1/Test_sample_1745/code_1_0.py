import sympy as sp

def solve_edl_potential():
    """
    This function derives and prints the expression for the Electrical Double-Layer (EDL)
    potential distribution in a parallel-plate microchannel based on the given problem.
    """
    # 1. Define the symbolic variables
    # y: spatial coordinate across the channel height
    # z1: base zeta potential at the bottom surface (j=1)
    # beta: slip length
    # k: Debye-Huckel parameter
    # H: height of the microchannel
    y, z1, beta, k, H = sp.symbols('y z_1 beta k H', real=True, positive=True)
    
    # psi(y) is the EDL potential function we want to find
    psi = sp.Function('psi')(y)

    # 2. Define the slip-dependent zeta potential at the bottom surface (z_a1)
    # According to the problem, z_aj = z_j(1 + beta*k).
    # Since z_2 = 0, the potential at the top wall is 0.
    # The potential at the bottom wall is z_a1.
    z_a1 = z1 * (1 + beta * k)

    # 3. Construct the final expression based on the analytical solution
    # The solution to the linearized Poisson-Boltzmann equation with the given
    # boundary conditions psi(H/2) = 0 and psi(-H/2) = z_a1 is:
    # psi(y) = z_a1 * sinh(k*(H/2 - y)) / sinh(k*H)

    # Numerator of the solution: sinh(k * (H/2 - y))
    numerator = sp.sinh(k * (sp.Rational(1, 2) * H - y))

    # Denominator of the solution: sinh(k * H)
    denominator = sp.sinh(k * H)

    # The full expression for the EDL potential distribution
    solution_expression = z_a1 * numerator / denominator

    # 4. Create and print the final equation
    final_equation = sp.Eq(psi, solution_expression)

    print("The derived expression for the Electrical double-layer potential distribution psi(y) is:")
    sp.pprint(final_equation, use_unicode=True)

    # As requested, outputting the numbers present in the final equation.
    print("\nThe numbers in the final equation are 1 and 2 (from the terms '1 + ...' and 'H/2').")


solve_edl_potential()

<<<psi(y) = z_1*(1 + beta*k)*sinh(k*(H/2 - y))/sinh(k*H)>>>