import sympy

def solve_linearization():
    """
    Solves for the linearization of the given predator-prey system
    at its non-trivial equilibrium point.
    """
    # Define symbolic variables
    S, F = sympy.symbols('S F')
    a, b, c, h, m = sympy.symbols('a b c h m')

    # Define the system of differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Substitute the given parameter values
    params = {a: 2, b: 1, c: 1, h: 1, m: 1}
    dS_dt_sub = dS_dt.subs(params)
    dF_dt_sub = dF_dt.subs(params)

    # Find the non-trivial equilibrium point (S > 0, F > 0)
    # We solve the system dS/dt = 0, dF/dt = 0
    # The equations derived in the explanation are S - F = 0 and 2 - F - S = 0
    equilibrium_sols = sympy.solve([S - F, 2 - F - S], (S, F))
    Se, Fe = equilibrium_sols[S], equilibrium_sols[F]
    
    # Create the Jacobian matrix
    system_matrix = sympy.Matrix([dS_dt_sub, dF_dt_sub])
    jacobian_matrix = system_matrix.jacobian([S, F])
    
    # Evaluate the Jacobian at the equilibrium point
    A_matrix = jacobian_matrix.subs({S: Se, F: Fe})
    
    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]
    
    # The constant terms in linearization around an equilibrium are zero.
    b11 = 0
    b22 = 0
    
    # Print the coefficients
    print(f"The non-trivial equilibrium point is (S, F) = ({Se}, {Fe})")
    print("\nThe linearized system is x' = Ax + b, where A is the Jacobian at the equilibrium point and b is the zero vector.")
    print("The coefficients are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_linearization()