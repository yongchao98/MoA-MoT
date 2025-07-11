import sympy

def solve_linearization():
    """
    Solves for the linearization of the Leslie's predator-prey system
    at the non-trivial equilibrium point.
    """
    # 1. Define variables and parameters
    S, F = sympy.symbols('S F')
    
    # Given parameters
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of ODEs
    # dS/dt = S * (h - m * S / F)
    # dF/dt = F * (a - b * F - c * S)
    f_S = S * (h - m * S / F)
    g_F = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium solution (Se, Fe)
    # For non-trivial solutions, S > 0 and F > 0.
    # From dS/dt = 0, we get h - m*S/F = 0 => S = (h/m)*F
    # Substitute S into dF/dt = 0: a - b*F - c*(h/m)*F = 0
    # F * (a - (b + c*h/m)) = 0. Since F > 0:
    # a = (b + c*h/m) * F
    Fe = a / (b + c * h / m)
    Se = (h / m) * Fe
    
    equilibrium_point = {S: Se, F: Fe}

    # 4. Compute the Jacobian matrix of the system
    system_matrix = sympy.Matrix([f_S, g_F])
    variables_matrix = sympy.Matrix([S, F])
    jacobian_matrix = system_matrix.jacobian(variables_matrix)
    
    # 5. Evaluate the Jacobian at the equilibrium point to find a_ij
    A_matrix = jacobian_matrix.subs(equilibrium_point)
    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # 6. The b_ij terms are the functions evaluated at the equilibrium, which are 0.
    b11 = f_S.subs(equilibrium_point)
    b22 = g_F.subs(equilibrium_point)

    # 7. Print the results.
    # The final equation format is:
    # x'(t) = a11*x(t) + a12*y(t) + b11
    # y'(t) = a21*x(t) + a22*y(t) + b22
    print(f"The non-trivial equilibrium point (S, F) is ({Se}, {Fe}).")
    print("The coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    solve_linearization()
