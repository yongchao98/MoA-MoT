import sympy as sp

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define the variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Step 2: Define the system of ODEs
    # dS/dt = S*(h - m*S/F)
    # dF/dt = F*(a - b*F - c*S)
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Step 3: Find the non-trivial equilibrium point (S0, F0)
    # We solve dS/dt = 0 and dF/dt = 0 for S > 0, F > 0.
    # h - m*S/F = 0  => 1 - S/F = 0 => S = F
    # a - b*F - c*S = 0 => 2 - F - S = 0
    # Substituting S = F into the second equation: 2 - F - F = 0 => 2*F = 2 => F = 1
    # Thus, S = 1.
    S0, F0 = 1, 1
    equilibrium_point = {S: S0, F: F0}

    # Step 4: Compute the Jacobian matrix of the system
    system_matrix = sp.Matrix([dS_dt, dF_dt])
    variable_matrix = sp.Matrix([S, F])
    jacobian_matrix = system_matrix.jacobian(variable_matrix)

    # Step 5: Evaluate the Jacobian at the equilibrium point to find a_ij
    jacobian_at_equilibrium = jacobian_matrix.subs(equilibrium_point)

    a11 = jacobian_at_equilibrium[0, 0]
    a12 = jacobian_at_equilibrium[0, 1]
    a21 = jacobian_at_equilibrium[1, 0]
    a22 = jacobian_at_equilibrium[1, 1]

    # Step 6: Determine the constant terms b_ij
    # The linearization around an equilibrium point results in a homogeneous system,
    # as the function value at the equilibrium point is zero.
    # Thus, the constant terms are zero.
    b11 = 0
    b22 = 0

    # Print the results
    print(f"The coefficients for the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    solve_linearization()