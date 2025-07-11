import sympy

def solve_linearization():
    """
    This function performs the linearization of the Leslie's predator-prey system
    at its non-trivial equilibrium point and prints the coefficients.
    """
    # Define symbolic variables
    S, F = sympy.symbols('S F')

    # Assign parameter values
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    # dS/dt = G(S, F)
    # dF/dt = H(S, F)
    G = S * (h - m * S / F)
    H = F * (a - b * F - c * S)

    # Solve for the non-trivial equilibrium point (S > 0, F > 0)
    # The equations are G = 0 and H = 0
    # For S>0, F>0, these simplify to:
    # h - m*S/F = 0  => 1 - S/F = 0 => S = F
    # a - b*F - c*S = 0 => 2 - F - S = 0
    # Substituting S=F into the second equation: 2 - F - F = 0 => 2*F = 2 => F = 1
    # Thus, S = 1.
    Se, Fe = 1, 1

    # Compute the Jacobian matrix of the system [G, H] with respect to [S, F]
    system_vector = sympy.Matrix([G, H])
    variables_vector = sympy.Matrix([S, F])
    J = system_vector.jacobian(variables_vector)

    # Substitute the equilibrium point into the Jacobian matrix to find the coefficient matrix A
    A = J.subs([(S, Se), (F, Fe)])

    # Extract the coefficients a_ij from the matrix A
    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    # For linearization at an equilibrium point, the constant vector B is zero.
    b11 = 0
    b22 = 0

    # Print the coefficients
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

solve_linearization()
<<<a_11=-1, a_12=1, a_21=-1, a_22=-1, b_11=0, b_22=0>>>