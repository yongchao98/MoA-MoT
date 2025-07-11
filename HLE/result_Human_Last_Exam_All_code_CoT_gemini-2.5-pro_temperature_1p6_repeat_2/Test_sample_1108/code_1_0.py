import sympy

def solve_linearization():
    """
    Solves for the linearization of the Leslie's predator-prey system at the non-trivial equilibrium.
    """
    # 1. Define symbolic variables and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium point (S > 0, F > 0)
    # We solve the system dS_dt = 0 and dF_dt = 0
    # From dS_dt = 0 (since S!=0), we have h - m*S/F = 0 -> S = (h/m)*F
    # Substitute S into dF_dt = 0 (since F!=0), we have a - b*F - c*(h/m)*F = 0
    # F * (b + c*h/m) = a -> F = a / (b + c*h/m)
    Fe_val = a / (b + c * h / m)
    Se_val = (h / m) * Fe_val
    
    equilibrium_point = {S: Se_val, F: Fe_val}

    # 4. Compute the Jacobian matrix of the system
    # The matrix entries are the coefficients a_ij
    J = sympy.Matrix([[sympy.diff(dS_dt, S), sympy.diff(dS_dt, F)],
                      [sympy.diff(dF_dt, S), sympy.diff(dF_dt, F)]])

    # 5. Evaluate the Jacobian at the equilibrium point
    A_matrix = J.subs(equilibrium_point)

    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # 6. The linearization at an equilibrium point is a homogeneous system, so b_ij are 0.
    b11 = 0
    b22 = 0

    # 7. Print the results
    print(f"The non-trivial equilibrium point (S_e, F_e) is ({Se_val}, {Fe_val})")
    print("\nThe linearized system is:")
    print("x'(t) = (a11 * x(t)) + (a12 * y(t)) + b11")
    print("y'(t) = (a21 * x(t)) + (a22 * y(t)) + b22")
    print("\nThe coefficients are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

solve_linearization()
<<<a_11 = -1, a_12 = 1, a_21 = -1, a_22 = -1, b_11 = 0, b_22 = 0>>>