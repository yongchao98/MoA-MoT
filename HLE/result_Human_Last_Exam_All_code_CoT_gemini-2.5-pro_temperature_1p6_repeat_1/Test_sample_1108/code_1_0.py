import sympy as sp

def solve_linearization():
    """
    This function solves for the linearization of the given predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define symbolic variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of differential equations
    # S' = f(S, F)
    # F' = g(S, F)
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium point (Se, Fe)
    # This is where dS/dt = 0 and dF/dt = 0, for S > 0, F > 0.
    # This simplifies the equations to:
    # h - m*S/F = 0
    # a - b*F - c*S = 0
    eq1 = sp.Eq(h - m * S / F, 0)
    eq2 = sp.Eq(a - b * F - c * S, 0)

    # Solve the system for S and F
    solution = sp.solve((eq1, eq2), (S, F))
    Se = solution[S]
    Fe = solution[F]

    # 4. Compute the Jacobian matrix of the system
    # J = [[df/dS, df/dF], [dg/dS, dg/dF]]
    system_matrix = sp.Matrix([dS_dt, dF_dt])
    variables_matrix = sp.Matrix([S, F])
    J = system_matrix.jacobian(variables_matrix)

    # 5. Evaluate the Jacobian at the equilibrium point to find the A matrix
    # A = [[a11, a12], [a21, a22]]
    A_matrix = J.subs([(S, Se), (F, Fe)])

    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # 6. Determine the constant vector b
    # This is the value of the system at the equilibrium point, which is [0, 0].
    b11 = dS_dt.subs([(S, Se), (F, Fe)])
    b22 = dF_dt.subs([(S, Se), (F, Fe)])

    # 7. Print the results
    print(f"The coefficients for the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == "__main__":
    solve_linearization()
