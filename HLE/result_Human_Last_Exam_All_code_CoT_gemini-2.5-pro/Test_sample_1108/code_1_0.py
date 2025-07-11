import sympy as sp

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at the non-trivial equilibrium point.
    """
    # Define symbolic variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = (2, 1, 1, 1, 1)

    # Define the system of differential equations
    # dS/dt = g(S, F)
    # dF/dt = f(S, F)
    g = S * (h - m * S / F)
    f = F * (a - b * F - c * S)

    # Find the non-trivial equilibrium point (S0, F0)
    # We solve the system g=0 and f=0 for S>0, F>0
    # From g=0 => h - m*S/F = 0 => S = (h/m)*F
    # Substitute into f=0 => a - b*F - c*(h/m)*F = 0 => F*(b + c*h/m) = a
    # F0 = a / (b + c*h/m)
    # S0 = (h/m)*F0
    F0 = a / (b + c * h / m)
    S0 = (h / m) * F0
    
    eq_point = {S: S0, F: F0}

    # Compute the Jacobian matrix of the system [g, f] with respect to [S, F]
    J = sp.Matrix([g, f]).jacobian([S, F])
    
    # Evaluate the Jacobian matrix at the equilibrium point
    J_eq = J.subs(eq_point)

    # The coefficients a_ij are the elements of the evaluated Jacobian matrix
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]

    # In the linearization around an equilibrium point, the constant terms are zero.
    b11 = 0
    b22 = 0

    # Print the results
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_linearization()