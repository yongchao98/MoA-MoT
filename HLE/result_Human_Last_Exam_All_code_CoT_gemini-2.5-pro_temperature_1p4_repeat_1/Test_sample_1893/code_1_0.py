import sympy
from sympy import cos, sin, pi

def solve_neutralino_mass():
    """
    Computes the eigenvalue of the neutralino mass matrix under specific conditions.
    
    The problem describes a "dynamic enhancement" scenario where the photino and a Higgsino
    are pure mass eigenstates. This implies two conditions on the parameters of the matrix:
    1. M_1 = M_2
    2. beta = pi / 4

    Under these conditions, the four eigenvalues generally depend on M_1 and mu.
    The question asks for an eigenvalue that is NOT proportional to M_1, M_2, or mu.
    This points to a contradiction unless a special case is considered. We will investigate
    the limiting case where the adjustable mass parameters M_1, M_2, and mu are set to zero.
    """
    
    print("Analyzing the neutralino mass matrix under dynamic enhancement conditions.")
    print("This requires M_1 = M_2 and beta = pi/4.")
    print("To find an eigenvalue not dependent on M_1, M_2, or mu, we consider the special case where M_1 = M_2 = mu = 0.")
    
    # Define symbols
    M_Z = sympy.Symbol('M_Z')
    theta_W = sympy.Symbol('theta_W')
    
    # Set parameter values for the special case
    M1 = 0
    M2 = 0
    mu = 0
    beta = pi/4
    
    # Define matrix elements based on the parameters
    c = cos(theta_W)
    s = sin(theta_W)
    
    m11 = M1 * c**2 + M2 * s**2
    m12 = (M2 - M1) * s * c
    m22 = M1 * s**2 + M2 * c**2
    m33 = mu * sin(2 * beta)
    m34 = -mu * cos(2 * beta)
    m44 = -mu * sin(2 * beta)
    
    # Construct the simplified neutralino mass matrix
    M_N_simplified = sympy.Matrix([
        [m11, m12,   0,   0],
        [m12, m22, M_Z,   0],
        [  0, M_Z, m33, m34],
        [  0,   0, m34, m44]
    ])
    
    print("\nThe simplified mass matrix with M_1 = M_2 = mu = 0 is:")
    sympy.pprint(M_N_simplified)
    
    # Calculate and print the eigenvalues of the simplified matrix
    eigenvalues = M_N_simplified.eigenvals()
    print("\nThe eigenvalues of this matrix are:")
    # The eigenvals() method returns a dictionary of {eigenvalue: multiplicity}
    eigenvalue_list = list(eigenvalues.keys())
    sympy.pprint(eigenvalues)

    # Identify the required eigenvalue
    result_eigenvalue = None
    for val in eigenvalue_list:
        # Find the eigenvalue that is not zero (as M1, M2, mu are zero)
        if val != 0:
            result_eigenvalue = val
            break
            
    # The physical mass is the absolute value of the eigenvalue, but the question asks for the eigenvalue itself.
    # We select one of the non-zero solutions, M_Z or -M_Z. Let's state the positive value.
    final_answer = M_Z
    
    print(f"\nTwo of the eigenvalues are 0. The other two are {M_Z} and -{M_Z}.")
    print(f"The eigenvalue that is not proportional to M_1, M_2, or mu is thus {final_answer} (or -{M_Z}).")
    print("\nThe final equation is:")
    print(f"Eigenvalue = {final_answer}")


solve_neutralino_mass()
<<<M_Z>>>