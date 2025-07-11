import sympy as sp

def solve_linearization():
    """
    This function solves for the linearization of the given Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define symbolic variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of differential equations
    # dS/dt = G1(S, F)
    # dF/dt = G2(S, F)
    G1 = S * (h - m * S / F)
    G2 = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium point (S0, F0)
    # We solve the system G1=0, G2=0 for S > 0, F > 0.
    # This is equivalent to solving the non-linear parts equal to zero.
    eq1 = h - m * S / F
    eq2 = a - b * F - c * S
    
    # Use sympy's solver
    solution = sp.solve([eq1, eq2], (S, F), dict=True)
    
    # The non-trivial solution requires S > 0 and F > 0.
    # The solver returns a list of solutions; we take the first one.
    non_trivial_solution = solution[0]
    S0 = non_trivial_solution[S]
    F0 = non_trivial_solution[F]

    # 4. Compute the Jacobian matrix of the system
    G = sp.Matrix([G1, G2])
    variables = sp.Matrix([S, F])
    J = G.jacobian(variables)

    # 5. Evaluate the Jacobian at the equilibrium point (S0, F0)
    # This gives the matrix A in the linearized system x' = Ax + b
    A = J.subs({S: S0, F: F0})

    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    # 6. Determine the constant vector b
    # For linearization at an equilibrium point, the constant vector is zero.
    b11 = 0
    b22 = 0
    
    # Print the results
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    solve_linearization()