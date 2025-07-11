import sympy

def solve_linearization():
    """
    Solves for the linearization of the Leslie's predator-prey system
    at the non-trivial equilibrium point.
    """
    # Step 1: Define symbols and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point
    # We are looking for S > 0 and F > 0
    # From dS/dt = 0, since S != 0, we have h - m*S/F = 0 => S = (h/m)*F
    # From dF/dt = 0, since F != 0, we have a - b*F - c*S = 0
    # Substitute S into the second equation: a - b*F - c*(h/m)*F = 0
    # F * (b + c*h/m) = a => F0 = a / (b + c*h/m)
    F0 = a / (b + c * h / m)
    S0 = (h / m) * F0
    
    # Ensure S0 and F0 are numeric values for jacobian substitution
    S0 = float(S0)
    F0 = float(F0)

    # Step 3: Compute the Jacobian matrix
    f = dS_dt
    g = dF_dt
    
    J = sympy.Matrix([f, g]).jacobian([S, F])
    
    # Step 4: Evaluate the Jacobian at the equilibrium point (S0, F0)
    J_eq = J.subs({S: S0, F: F0})

    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # Step 5: Determine the constant terms
    # For linearization around an equilibrium point, the constant terms are zero.
    b11 = 0
    b22 = 0

    # Print the results
    print(f"The non-trivial equilibrium point is (S0, F0) = ({S0}, {F0})")
    print("\nThe coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    solve_linearization()

<<<-1, 1, -1, -1, 0, 0>>>