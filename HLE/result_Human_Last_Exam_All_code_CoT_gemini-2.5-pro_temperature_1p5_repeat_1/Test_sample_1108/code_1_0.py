import sympy as sp

def find_linearization_coeffs():
    """
    This function calculates the coefficients of the linearized system for the given
    Leslie's predator-prey model at its non-trivial equilibrium point.
    """
    # 1. Define symbolic variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = (2, 1, 1, 1, 1)

    # Define the system of differential equations with given parameters
    # dS/dt = S(h - m*S/F)
    # dF/dt = F(a - b*F - c*S)
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # 2. Find the non-trivial equilibrium point (S > 0, F > 0)
    # The equations to solve are:
    # h - m*S/F = 0  => 1 - S/F = 0 => S = F
    # a - b*F - c*S = 0 => 2 - F - S = 0
    # Substituting S = F into the second equation: 2 - F - F = 0 => 2*F = 2 => F = 1.
    # Therefore, S = 1.
    equilibrium_point = {S: 1, F: 1}
    S_e = equilibrium_point[S]
    F_e = equilibrium_point[F]

    # 3. Compute the Jacobian matrix
    system_vector = sp.Matrix([dS_dt, dF_dt])
    variables_vector = sp.Matrix([S, F])
    jacobian_matrix = system_vector.jacobian(variables_vector)

    # 4. Evaluate the Jacobian at the equilibrium point
    # This gives the matrix A = [[a11, a12], [a21, a22]]
    A = jacobian_matrix.subs(equilibrium_point)

    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    # 5. The constant terms b11 and b22 are 0 for linearization at an equilibrium point.
    b11 = 0
    b22 = 0

    # Print the coefficients
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    find_linearization_coeffs()
<<<a_11 = -1, a_12 = 1, a_21 = -1, a_22 = -1, b_11 = 0, b_22 = 0>>>