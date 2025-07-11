import math

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define parameters
    a = 2.0
    b = 1.0
    c = 1.0
    h = 1.0
    m = 1.0

    # Step 2: Find the non-trivial equilibrium point (Se, Fe)
    # The system is dS/dt = 0 and dF/dt = 0.
    # For S > 0, F > 0:
    # From dS/dt = S(h - m*S/F) = 0, we get h = m*S/F => S = (h/m)*F.
    # From dF/dt = F(a - b*F - c*S) = 0, we get a = b*F + c*S.
    # Substitute S into the second equation: a = b*F + c*(h/m)*F = F*(b + c*h/m).
    Fe = a / (b + c * h / m)
    Se = (h / m) * Fe

    # Step 3: Compute the coefficients a_ij by evaluating the Jacobian matrix at (Se, Fe).
    # The partial derivatives are:
    # d/dS (S(h - mS/F)) = h - 2*m*S/F
    # d/dF (S(h - mS/F)) = m*S^2/F^2
    # d/dS (F(a - bF - cS)) = -c*F
    # d/dF (F(a - bF - cS)) = a - 2*b*F - c*S

    a11 = h - 2 * m * Se / Fe
    a12 = m * (Se**2) / (Fe**2)
    a21 = -c * Fe
    a22 = a - 2 * b * Fe - c * Se

    # Step 4: The constant terms b_ij are zero for a linearization at an equilibrium point.
    b11 = 0.0
    b22 = 0.0
    
    # Step 5: Print the results.
    print("The coefficients of the linearized system are:")
    # Using math.isclose to handle potential floating point inaccuracies for integer conversion
    print(f"a_11 = {int(a11) if math.isclose(a11, round(a11)) else a11}")
    print(f"a_12 = {int(a12) if math.isclose(a12, round(a12)) else a12}")
    print(f"a_21 = {int(a21) if math.isclose(a21, round(a21)) else a21}")
    print(f"a_22 = {int(a22) if math.isclose(a22, round(a22)) else a22}")
    print(f"b_11 = {int(b11)}")
    print(f"b_22 = {int(b22)}")

solve_linearization()