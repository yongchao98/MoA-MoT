import sympy

def solve_edl_potential():
    """
    This function programmatically constructs and prints the expression for the
    Electrical Double-Layer (EDL) potential distribution based on the problem's conditions.
    """
    # Define the symbolic variables from the problem description
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H')

    # The slip-dependant zeta potential at the bottom surface (y=0) is z_a1.
    # The potential at the top surface (y=H) is 0 since z_2=0.
    z_a1 = z_1 * (1 + beta * k)

    # The final expression for the EDL potential psi(y) is derived from the
    # linearized Poisson-Boltzmann equation with the given boundary conditions.
    # As derived in the plan, the solution is:
    # psi(y) = z_a1 * sinh(k*(H-y)) / sinh(kH)
    final_expression = z_a1 * sympy.sinh(k * (H - y)) / sympy.sinh(k * H)

    # We print the derived expression. Sympy helps in formatting it correctly.
    # Note: Sympy might rearrange terms, e.g., (1 + beta*k) to (beta*k + 1)
    # and (H - y) to (-y + H), which is mathematically equivalent.
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # Create the full equation psi(y) = expression for a clear output
    psi_y = sympy.Function('psi')(y)
    full_equation = sympy.Eq(psi_y, final_expression)
    print(str(full_equation))
    
    # Per the instructions, return the final answer in the specified format.
    # The "answer" is the expression for psi(y).
    print(f"\n<<<{str(final_expression)}>>>")

solve_edl_potential()