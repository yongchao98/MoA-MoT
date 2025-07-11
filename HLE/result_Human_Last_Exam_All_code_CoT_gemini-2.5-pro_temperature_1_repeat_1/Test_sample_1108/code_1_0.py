def solve_leslie_linearization():
    """
    Solves for the linearization of the Leslie's predator-prey system
    at the non-trivial equilibrium point.
    """
    # Step 1: Define parameters
    a = 2.0
    b = 1.0
    c = 1.0
    h = 1.0
    m = 1.0

    # Step 2: Find the non-trivial equilibrium point (Se, Fe)
    # The equations at equilibrium are:
    # h - m*S/F = 0  => S = (h/m)*F
    # a - b*F - c*S = 0
    # Substituting S into the second equation: a - b*F - c*(h/m)*F = 0
    # a = F * (b + c*h/m) => F = a / (b + c*h/m)
    Fe = a / (b + c * h / m)
    Se = (h / m) * Fe

    # Step 3 & 4: Define the partial derivatives (Jacobian elements)
    # and evaluate them at the equilibrium point (Se, Fe)
    # a11 = d/dS (hS - m*S^2/F) = h - 2*m*S/F
    a11 = h - 2 * m * Se / Fe
    
    # a12 = d/dF (hS - m*S^2/F) = m*S^2/F^2
    a12 = m * (Se**2) / (Fe**2)
    
    # a21 = d/dS (aF - b*F^2 - c*F*S) = -c*F
    a21 = -c * Fe
    
    # a22 = d/dF (aF - b*F^2 - c*F*S) = a - 2*b*F - c*S
    a22 = a - 2 * b * Fe - c * Se

    # Step 5: The constant terms in the linearization around an equilibrium are zero.
    b11 = 0.0
    b22 = 0.0

    # Print the results
    print(f"The non-trivial equilibrium point is (S_e, F_e) = ({Se}, {Fe})")
    print("\nThe coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

    # Print the final equation with the numbers
    print("\nThe linearized system in the form x' = Ax + b is:")
    print(f"x'(t) = ({a11})*x(t) + ({a12})*y(t) + ({b11})")
    print(f"y'(t) = ({a21})*x(t) + ({a22})*y(t) + ({b22})")

# Execute the function
solve_leslie_linearization()
