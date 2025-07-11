import math

def solve_linearization():
    """
    This function solves for the linearization of the Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. System parameters
    a = 2.0
    b = 1.0
    c = 1.0
    h = 1.0
    m = 1.0

    # 2. Find the non-trivial equilibrium point (Se, Fe)
    # From dS/dt = S(h - m*S/F) = 0 and S > 0, we have h - m*S/F = 0 => S/F = h/m.
    # From dF/dt = F(a - b*F - c*S) = 0 and F > 0, we have a - b*F - c*S = 0.
    # Substitute S = (h/m)*F into the second equation: a - b*F - c*(h/m)*F = 0.
    # This gives Fe = a / (b + c*h/m) and Se = (h/m)*Fe.
    Fe = a / (b + c * h / m)
    Se = (h / m) * Fe

    # 3. & 4. Calculate Jacobian elements at the equilibrium point
    # The system functions are:
    # f(S, F) = S*(h - m*S/F)
    # g(S, F) = F*(a - b*F - c*S)
    # The partial derivatives are:
    # df/dS = h - 2*m*S/F
    # df/dF = m*S^2/F^2
    # dg/dS = -c*F
    # dg/dF = a - 2*b*F - c*S
    
    # Evaluate at the equilibrium point (Se, Fe)
    a11 = h - 2 * m * Se / Fe
    a12 = m * (Se**2) / (Fe**2)
    a21 = -c * Fe
    a22 = a - 2 * b * Fe - c * Se

    # 5. The b vector components are zero at an equilibrium point.
    b11 = 0.0
    b22 = 0.0
    
    # Print the values as requested
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

solve_linearization()