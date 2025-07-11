import sympy

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define the system and parameters symbolically
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Right-hand side functions of the differential equations
    f_SF = S * (h - m * S / F)
    g_SF = F * (a - b * F - c * S)

    # 2. Find the non-trivial equilibrium point (S_e, F_e)
    # We solve f_SF=0 and g_SF=0 for S>0, F>0
    # From f_SF=0 (with S>0), we get h - m*S/F = 0 => F = (m/h)*S => F = S
    # Substitute F=S into g_SF=0 (with F>0): a - b*S - c*S = 0 => a = (b+c)*S
    S_e = a / (b + c)
    F_e = S_e
    
    # 3. Calculate the Jacobian matrix
    J = sympy.Matrix([f_SF, g_SF]).jacobian([S, F])
    
    # 4. Evaluate the Jacobian at the equilibrium point to get a_ij
    J_eq = J.subs({S: S_e, F: F_e})
    
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]
    
    # 5. Determine b_ij, which are the function values at equilibrium (must be 0)
    b11 = f_SF.subs({S: S_e, F: F_e})
    b22 = g_SF.subs({S: S_e, F: F_e})

    # Print the results
    print(f"The non-trivial equilibrium point is (S_e, F_e) = ({S_e}, {F_e})")
    print("\nThe coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

solve_linearization()