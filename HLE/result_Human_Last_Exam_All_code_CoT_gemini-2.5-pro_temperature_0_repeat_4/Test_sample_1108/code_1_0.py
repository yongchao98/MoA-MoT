import sympy

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define the system with given parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Right-hand side functions of the differential equations
    f = S * (h - m * S / F)
    g = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point (S_e, F_e)
    # We solve f=0 and g=0 for S>0, F>0.
    # sympy.solve can find the solutions to the system of equations.
    equilibria = sympy.solve([f, g], [S, F], dict=True)
    
    # Filter for the non-trivial solution where S > 0 and F > 0
    non_trivial_eq = None
    for eq in equilibria:
        if eq[S] > 0 and eq[F] > 0:
            non_trivial_eq = eq
            break
    
    if non_trivial_eq is None:
        print("Error: Could not find a non-trivial equilibrium point.")
        return

    S_e = non_trivial_eq[S]
    F_e = non_trivial_eq[F]

    # Step 3: Compute the Jacobian matrix
    # The Jacobian matrix contains the partial derivatives of f and g with respect to S and F.
    J = sympy.Matrix([f, g]).jacobian([S, F])

    # Step 4: Evaluate the Jacobian at the equilibrium point
    # This gives the matrix of coefficients for the linearized system.
    J_eq = J.subs([(S, S_e), (F, F_e)])
    
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # Step 5: Determine the constant terms
    # For a linearization around an equilibrium point, the constant vector is zero
    # because the functions f and g are zero at this point.
    b11 = 0
    b22 = 0

    # Step 6: Print the results as requested
    print(f"The non-trivial equilibrium point is (S_e, F_e) = ({S_e}, {F_e}).")
    print("The coefficients for the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")
    
    print("\nThe final linearized system is:")
    print(f"x'(t) = ({a11})*x(t) + ({a12})*y(t) + ({b11})")
    print(f"y'(t) = ({a21})*x(t) + ({a22})*y(t) + ({b22})")

if __name__ == '__main__':
    solve_linearization()