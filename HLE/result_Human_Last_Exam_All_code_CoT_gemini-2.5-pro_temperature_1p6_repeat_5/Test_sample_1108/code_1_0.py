import sympy

def solve_linearization():
    """
    This function finds the linearization of the Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define the system with given parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1
    
    # System equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point (S > 0, F > 0)
    # Solve dS/dt = 0 and dF/dt = 0
    # From dS/dt = S * (1 - S/F) = 0, for S > 0, we have 1 - S/F = 0 => S = F
    # Substitute S = F into dF/dt = 0:
    # F * (2 - F - F) = 0 => F * (2 - 2F) = 0
    # For F > 0, we have 2 - 2F = 0 => F = 1
    # Thus, the equilibrium point is S_0 = 1, F_0 = 1
    S0, F0 = 1, 1

    # Step 3: Compute the Jacobian matrix
    f = dS_dt
    g = dF_dt
    
    J = sympy.Matrix([
        [sympy.diff(f, S), sympy.diff(f, F)],
        [sympy.diff(g, S), sympy.diff(g, F)]
    ])

    # Step 4: Evaluate the Jacobian at the equilibrium point (S0, F0)
    J_eq = J.subs([(S, S0), (F, F0)])
    
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # Step 5: Define the constant terms
    # For a linearization at an equilibrium point, the constant terms are zero.
    b11 = 0
    b22 = 0

    # Print the results
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

solve_linearization()
<<<[-1, 1, -1, -1, 0, 0]>>>