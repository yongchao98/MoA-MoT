import sympy as sp

def solve_linearization():
    """
    This function calculates the linearization of the given predator-prey system
    at its non-trivial equilibrium point and prints the coefficients.
    """
    # Step 1: Define the system variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # The system of differential equations
    # dS/dt = S(h - m*S/F)
    # dF/dt = F(a - b*F - c*S)
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point (S0, F0)
    # We solve dS/dt = 0 and dF/dt = 0 for S > 0 and F > 0.
    # The non-trivial solution is found by solving the system:
    # h - m*S/F = 0  => 1 - S/F = 0 => S = F
    # a - b*F - c*S = 0 => 2 - F - S = 0
    # Substituting S = F into the second equation gives 2 - 2F = 0, so F = 1.
    # Thus, S = 1.
    S0, F0 = 1, 1

    # Step 3: Compute the Jacobian matrix of the system
    # The Jacobian contains the partial derivatives of the system's functions.
    system_matrix = sp.Matrix([dS_dt, dF_dt])
    variables_matrix = sp.Matrix([S, F])
    J = system_matrix.jacobian(variables_matrix)

    # Step 4: Evaluate the Jacobian at the equilibrium point (S0, F0)
    # This gives the matrix A in the linearized system.
    A_matrix = J.subs({S: S0, F: F0})

    # Extract the coefficients a_ij from the resulting matrix
    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # Step 5: Determine the constant terms b_ij
    # For a linearization at an equilibrium point, the system is homogeneous,
    # meaning the constant vector is zero.
    b11 = 0
    b22 = 0

    # Print the values of each coefficient
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_linearization()