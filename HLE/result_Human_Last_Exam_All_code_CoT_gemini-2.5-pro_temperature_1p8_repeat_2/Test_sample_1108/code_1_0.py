import sympy

def solve_linearization():
    """
    This function solves for the linearization of the given Leslie's predator-prey system.
    1. It defines the system of ODEs with the specified parameters.
    2. It finds the non-trivial equilibrium point.
    3. It computes the Jacobian matrix of the system.
    4. It evaluates the Jacobian at the equilibrium point to find the coefficients a_ij.
    5. It sets the constant terms b_ii, which are zero for linearization at an equilibrium.
    6. It prints the final coefficients.
    """
    # Define symbolic variables
    S, F = sympy.symbols('S F')
    
    # Parameters for the system
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    # dS/dt = f(S, F)
    f = S * (h - m * S / F)
    # dF/dt = g(S, F)
    g = F * (a - b*F - c*S)

    # Solve for the non-trivial equilibrium point (S > 0, F > 0)
    # The solver returns all equilibrium points, we filter for the one where S > 0 and F > 0.
    equilibria = sympy.solve([f, g], (S, F), dict=True)
    non_trivial_eq = None
    for sol in equilibria:
        if sol[S] > 0 and sol[F] > 0:
            non_trivial_eq = sol
            break
    
    if non_trivial_eq is None:
        print("Non-trivial equilibrium not found.")
        return

    S0 = non_trivial_eq[S]
    F0 = non_trivial_eq[F]

    # Compute the Jacobian matrix of the system [f, g] with respect to [S, F]
    system_matrix = sympy.Matrix([f, g])
    J = system_matrix.jacobian([S, F])

    # Evaluate the Jacobian at the equilibrium point (S0, F0)
    A_matrix = J.subs({S: S0, F: F0})
    
    a11 = A_matrix[0, 0]
    a12 = A_matrix[0, 1]
    a21 = A_matrix[1, 0]
    a22 = A_matrix[1, 1]

    # For linearization at an equilibrium point, the constant term is zero.
    b11 = 0
    b22 = 0

    # Print the resulting coefficients
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == "__main__":
    solve_linearization()