def calculate_linearization_coefficients():
    """
    Calculates the coefficients of the linearized Leslie's predator-prey system
    at the non-trivial equilibrium point.
    """
    # System parameters
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Step 1: Find the non-trivial equilibrium point (Se, Fe)
    # From dS/dt = 0, we get S = (h/m)*F. With given params, S = F.
    # From dF/dt = 0, we get a - b*F - c*S = 0.
    # Substituting S=F: a - b*F - c*F = 0 => a = (b+c)*F
    Fe = a / (b + c)
    Se = Fe
    
    # Step 2: Calculate the coefficients a_ij of the Jacobian matrix at (Se, Fe)
    # Partial derivatives:
    # df/dS = h - 2*m*S/F
    # df/dF = m*S^2/F^2
    # dg/dS = -c*F
    # dg/dF = a - 2*b*F - c*S
    
    a11 = h - 2 * m * Se / Fe
    a12 = m * Se**2 / Fe**2
    a21 = -c * Fe
    a22 = a - 2 * b * Fe - c * Se
    
    # Step 3: The constant terms b_ij are zero by definition of the equilibrium point.
    b11 = 0
    b22 = 0
    
    # Step 4: Print the results
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

calculate_linearization_coefficients()