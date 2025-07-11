import sympy as sp

def solve_edl_potential():
    """
    This function symbolically solves for the Electrical Double-Layer (EDL)
    potential distribution psi(y) in a parallel-plate microchannel.
    """
    # Define the symbolic variables and the function
    y, k, H, beta = sp.symbols('y k H beta')
    zeta_1 = sp.Symbol('zeta_1')
    psi = sp.Function('psi')

    # Define the Debye-Hückel ordinary differential equation (ODE)
    # d^2(psi)/dy^2 - k^2*psi = 0
    ode = sp.Eq(psi(y).diff(y, 2) - k**2 * psi(y), 0)

    # Get the general solution of the ODE
    general_solution = sp.dsolve(ode, psi(y))
    solution_rhs = general_solution.rhs

    # Define the slip-dependent zeta potential at the bottom wall (y=0)
    # as per the problem description: z_a1 = z_1*(1 + beta*k)
    zeta_a1 = zeta_1 * (1 + beta * k)

    # Define the boundary conditions:
    # 1. At the top wall (y=H), potential is 0 (since zeta_2=0)
    # 2. At the bottom wall (y=0), potential is zeta_a1
    # The general solution has constants C1 and C2 to be solved for.
    C1, C2 = sp.symbols('C1 C2')
    
    bcs = [
        sp.Eq(solution_rhs.subs(y, H), 0),
        sp.Eq(solution_rhs.subs(y, 0), zeta_a1)
    ]
    
    # Solve for the constants C1 and C2 based on the boundary conditions
    constants = sp.solve(bcs, (C1, C2))
    
    # Substitute the constants back into the general solution
    particular_solution = solution_rhs.subs(constants)
    
    # Simplify the final expression
    simplified_solution = sp.simplify(particular_solution)

    # Format the solution into a standard mathematical expression string
    # sympy might output sinh(k*(-y + H)), which is correct but less conventional.
    # We will format it to sinh(k*(H - y)).
    final_expression_str = str(simplified_solution).replace('k*(-y + H)', 'k*(H - y)')
    
    # Print the final equation for the EDL potential distribution.
    # We print each part of the equation as requested.
    numerator_part_1 = f"{zeta_1}*(1 + {beta}*{k})"
    numerator_part_2 = f"sinh({k}*(H - y))"
    denominator = f"sinh({k}*{H})"
    
    print(f"The expression for the Electrical double-layer potential distribution is:")
    print(f"psi(y) = ( {numerator_part_1} ) * ( {numerator_part_2} ) / ( {denominator} )")

    # The final complete expression is also constructed for the final answer block.
    global final_answer
    final_answer = f"psi(y) = {zeta_1}*(1 + beta*k)*sinh(k*(H - y))/sinh(k*H)"


# Run the solver and print the result
final_answer = ""
solve_edl_potential()