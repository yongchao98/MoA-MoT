import sympy

def solve_linearization():
    """
    This function calculates the linearization of the Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define symbolic variables and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of differential equations
    # dS/dt = f(S, F)
    # dF/dt = g(S, F)
    f_S_F = S * (h - m * S / F)
    g_S_F = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium point (S_eq, F_eq)
    # Solve dS/dt = 0 and dF/dt = 0 for S > 0, F > 0
    # sympy.solve might return multiple solutions, we need to find the one where S and F are positive.
    equilibria = sympy.solve((f_S_F, g_S_F), (S, F), dict=True)
    
    non_trivial_eq = None
    for sol in equilibria:
        if sol[S] > 0 and sol[F] > 0:
            non_trivial_eq = sol
            break
            
    if non_trivial_eq is None:
        print("Non-trivial equilibrium not found.")
        return

    S_eq = non_trivial_eq[S]
    F_eq = non_trivial_eq[F]

    # 4. Compute the Jacobian matrix of the system
    system_matrix = sympy.Matrix([f_S_F, g_S_F])
    variables_matrix = sympy.Matrix([S, F])
    J = system_matrix.jacobian(variables_matrix)

    # 5. Evaluate the Jacobian at the non-trivial equilibrium point
    J_eq = J.subs({S: S_eq, F: F_eq})

    # 6. Extract the coefficients a_ij for the linearized system
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # 7. For linearization around an equilibrium point, the constant terms are zero.
    b11 = 0
    b22 = 0

    # 8. Print the final equations and the coefficients
    print(f"The non-trivial equilibrium point (S, F) is ({S_eq}, {F_eq}).")
    print("The linearized system is given by:")
    print(f"x'(t) = ({a11})*x(t) + ({a12})*y(t) + ({b11})")
    print(f"y'(t) = ({a21})*x(t) + ({a22})*y(t) + ({b22})")
    print("\nThe values of the coefficients are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_linearization()