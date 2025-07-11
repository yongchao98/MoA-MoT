import sympy as sp

def solve_linearization():
    """
    This function computes the linearization of the given predator-prey system
    at its non-trivial equilibrium point.
    """
    # Define symbolic variables
    S, F = sp.symbols('S F')

    # Define system parameters
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system equations
    f_S = S * (h - m * S / F)
    g_F = F * (a - b * F - c * S)

    # Solve for the non-trivial equilibrium point (S > 0, F > 0)
    # The equations to solve come from the terms in parentheses
    equilibrium_eqs = [
        sp.Eq(h - m * S / F, 0),
        sp.Eq(a - b * F - c * S, 0)
    ]
    # We solve this simple system manually for clarity, though sympy.solve could be used.
    # From the first equation: S = F
    # Substituting into the second: 2 - S - S = 0  => 2*S = 2 => S = 1
    S0, F0 = 1, 1

    # Define the vector field for Jacobian calculation
    vector_field = sp.Matrix([f_S, g_F])

    # Compute the Jacobian matrix
    jacobian_matrix = vector_field.jacobian([S, F])

    # Substitute the equilibrium point into the Jacobian to get the matrix A
    A = jacobian_matrix.subs({S: S0, F: F0})

    # Extract the coefficients a_ij
    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    # The vector b is the value of the vector field at equilibrium, which is [0, 0]
    b11 = f_S.subs({S: S0, F: F0})
    b22 = g_F.subs({S: S0, F: F0})
    
    # Print the coefficients
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == "__main__":
    solve_linearization()