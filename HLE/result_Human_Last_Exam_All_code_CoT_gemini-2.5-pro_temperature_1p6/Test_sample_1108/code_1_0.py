import sympy

def find_linearization_coefficients():
    """
    This function performs the linearization of the given predator-prey system
    at its non-trivial equilibrium point and prints the coefficients.
    """
    # Define symbolic variables for the state variables and parameters
    S, F = sympy.symbols('S F', positive=True)
    a, b, c, h, m = sympy.symbols('a b c h m')

    # Substitute the given parameter values
    parameter_values = {a: 2, b: 1, c: 1, h: 1, m: 1}

    # Define the system of differential equations
    # dS/dt = f(S, F)
    # dF/dt = g(S, F)
    f_S_F = S * (h - m * S / F)
    g_S_F = F * (a - b * F - c * S)

    # Substitute the parameter values into the equations
    f_eq = f_S_F.subs(parameter_values)
    g_eq = g_S_F.subs(parameter_values)

    # Find the non-trivial equilibrium point (S_e, F_e) by solving
    # f(S, F) = 0 and g(S, F) = 0 for S > 0, F > 0.
    equilibrium_solutions = sympy.solve([f_eq, g_eq], (S, F), dict=True)

    # Filter for the non-trivial solution where S and F are positive
    non_trivial_equilibrium = None
    for sol in equilibrium_solutions:
        if sol[S] > 0 and sol[F] > 0:
            non_trivial_equilibrium = sol
            break
            
    if non_trivial_equilibrium is None:
        print("Could not find the non-trivial equilibrium point.")
        return

    S_e = non_trivial_equilibrium[S]
    F_e = non_trivial_equilibrium[F]

    # Compute the Jacobian matrix of the system
    system_matrix = sympy.Matrix([f_eq, g_eq])
    variables_matrix = sympy.Matrix([S, F])
    J = system_matrix.jacobian(variables_matrix)

    # Evaluate the Jacobian at the non-trivial equilibrium point
    J_at_equilibrium = J.subs({S: S_e, F: F_e})

    # Extract the coefficients a_ij from the Jacobian matrix
    a11 = J_at_equilibrium[0, 0]
    a12 = J_at_equilibrium[0, 1]
    a21 = J_at_equilibrium[1, 0]
    a22 = J_at_equilibrium[1, 1]

    # For linearization at an equilibrium point, the constant terms are zero.
    b11 = 0
    b22 = 0

    # Print the final coefficients as requested
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == "__main__":
    find_linearization_coefficients()