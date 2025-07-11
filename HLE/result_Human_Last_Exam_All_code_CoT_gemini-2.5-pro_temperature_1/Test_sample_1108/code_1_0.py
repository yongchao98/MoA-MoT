import sympy as sp

def solve_linearization():
    """
    This function solves for the linearization of the given Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Define symbolic variables
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    # f represents dS/dt, g represents dF/dt
    f = S * (h - m * S / F)
    g = F * (a - b * F - c * S)

    # Find the non-trivial equilibrium point (S > 0, F > 0)
    # We solve f = 0 and g = 0
    # The non-trivial solution comes from:
    # 1 - S/F = 0  => S = F
    # 2 - F - S = 0
    # Substituting S = F into the second equation gives 2 - 2F = 0, so F = 1.
    # Therefore, S = 1.
    S_eq, F_eq = 1, 1

    # Compute the Jacobian matrix of the system [f, g] with respect to [S, F]
    X = sp.Matrix([f, g])
    Y = sp.Matrix([S, F])
    J = X.jacobian(Y)

    # Evaluate the Jacobian at the equilibrium point (S_eq, F_eq)
    J_eq = J.subs({S: S_eq, F: F_eq})

    # The coefficients a_ij are the entries of the evaluated Jacobian matrix
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # The linearization at an equilibrium point results in a homogeneous system
    # of the form x' = J(x_eq) * x. Thus, the constant terms b_ij are zero.
    b11 = 0
    b22 = 0

    # Print the resulting coefficients
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

solve_linearization()
<<<a11 = -1, a12 = 1, a21 = -1, a22 = -1, b11 = 0, b22 = 0>>>