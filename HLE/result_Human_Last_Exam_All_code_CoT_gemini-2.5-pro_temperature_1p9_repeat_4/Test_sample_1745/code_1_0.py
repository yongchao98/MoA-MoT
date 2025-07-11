import sympy

def solve_edl_potential():
    """
    This function derives the expression for the Electrical Double-Layer (EDL)
    potential distribution in a parallel-plate microchannel based on the given
    problem description.
    """

    # 1. Define the symbolic variables representing the physical parameters.
    # y: spatial coordinate across the channel height
    # H: height of the microchannel
    # k: Debye-Huckel parameter
    # beta: slip length
    # z_1: zeta potential parameter for the bottom surface
    y, H, k, beta, z_1 = sympy.symbols('y H k beta z_1', real=True)
    psi = sympy.Function('psi')

    # 2. Define the linearized Poisson-Boltzmann equation: d^2(psi)/dy^2 - k^2*psi = 0
    diff_eq = sympy.Eq(psi(y).diff(y, 2) - k**2 * psi(y), 0)

    # 3. Define the boundary conditions based on the problem statement.
    # The coordinate system is centered, with walls at y = -H/2 and y = +H/2.
    # Bottom wall (y = -H/2): Potential is the slip-dependent zeta potential z_a1.
    zeta_a1 = z_1 * (1 + sympy.Symbol('beta') * sympy.Symbol('k'))
    
    # Top wall (y = +H/2): Potential is 0, since z_2 = 0.
    bcs = {psi(H/2): 0, psi(-H/2): zeta_a1}

    # 4. Solve the ODE with the specified boundary conditions using sympy's dsolve.
    solution = sympy.dsolve(diff_eq, ics=bcs)

    # The result is an Equation, so we extract the right-hand side.
    psi_expression = solution.rhs

    # 5. Simplify the expression for better readability.
    # The standard form for this solution uses hyperbolic sine functions.
    simplified_psi_expression = sympy.simplify(psi_expression)

    # 6. Print the final equation for the potential distribution.
    print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
    final_equation_str = f"{psi(y)} = {simplified_psi_expression}"
    print(final_equation_str)

    # Output the answer in the requested format.
    answer_expression = str(simplified_psi_expression)
    print(f"\n<<<{answer_expression}>>>")

if __name__ == '__main__':
    solve_edl_potential()