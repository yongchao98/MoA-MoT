import sympy

def solve_neutralino_mass():
    """
    Computes the special eigenvalue of the neutralino mass matrix.
    """
    # Define symbols for the parameters in the matrix
    M1, M2, mu, M_Z = sympy.symbols('M_1 M_2 mu M_Z', real=True)
    theta_W, beta = sympy.symbols('theta_W beta', real=True)

    # Define symbolic trigonometric functions for clarity
    cW = sympy.cos(theta_W)
    sW = sympy.sin(theta_W)
    cos2b = sympy.cos(2*beta)
    sin2b = sympy.sin(2*beta)

    # Define the neutralino mass matrix as given in the problem
    M_N = sympy.Matrix([
        [M1*cW**2 + M2*sW**2, (M2 - M1)*sW*cW, 0, 0],
        [(M2 - M1)*sW*cW, M1*sW**2 + M2*cW**2, M_Z, 0],
        [0, M_Z, mu*sin2b, -mu*cos2b],
        [0, 0, -mu*cos2b, -mu*sin2b]
    ])

    print("The neutralino mass matrix in the (PHOTINO, ZINO, HIGGSINO_A, HIGGSINO_B) basis is:")
    sympy.pprint(M_N)

    print("\nThe scenario of dynamic enhancement requires the photino and a higgsino state to decouple from the zino and the other higgsino.")
    print("This implies that the mixing terms between these two sectors must be zero.")
    print("The conditions for this are:")
    print("1. (M2 - M1) * sin(theta_W) * cos(theta_W) = 0  =>  M1 = M2")
    print("2. -mu * cos(2*beta) = 0                      =>  cos(2*beta) = 0  =>  beta = pi/4 (for 0 <= beta <= pi/2)")
    
    print("\nWe are asked for an eigenvalue that is not proportional to M1, M2, or mu.")
    print("This points to an eigenvalue that is a constant. We investigate a specific parameter choice that fulfills the enhancement conditions and the physical description.")
    print("Let's consider the case where the gaugino masses are zero: M1 = M2 = 0. This also fulfills M1 = M2.")
    print("With M1 = M2 = 0 and beta = pi/4, the matrix becomes:")

    # Substitute the special conditions into the matrix
    # M1 = M2 = 0 and beta = pi/4
    M_special = M_N.subs({
        M1: 0,
        M2: 0,
        beta: sympy.pi / 4
    })
    
    sympy.pprint(M_special)

    # The eigenvalues of M_special can be found. It is block diagonal.
    # The top-left element is 0, so one eigenvalue is 0.
    # The other 3x3 block is [[0, M_Z, 0], [M_Z, mu, 0], [0, 0, -mu]]
    # This 3x3 block has eigenvalues -mu, and the eigenvalues of [[0, M_Z], [M_Z, mu]].
    # The final eigenvalues are 0, -mu, and the solutions to lambda^2 - mu*lambda - M_Z^2 = 0.
    eigenvalues = M_special.eigenvals() # Returns a dict {eigenvalue: multiplicity}

    print("\nThe eigenvalues of this specific matrix are:")
    # Pretty print the list of eigenvalue expressions
    sympy.pprint(list(eigenvalues.keys()))
    
    # Identify the eigenvalue that is a numerical constant and does not depend on other parameters
    result_eigenvalue = None
    for val in eigenvalues.keys():
        # Check if the eigenvalue is a constant number
        if val.is_Number:
            result_eigenvalue = val
            break
            
    print("\nThe eigenvalue which is not proportional to M1, M2, or mu is the constant value:")
    print(result_eigenvalue)
    
    # We will print the final answer in the required format
    # In this case, the final equation is just the value of the eigenvalue.
    print("\nFinal Equation:")
    print("eigenvalue = ", end="")
    sympy.pprint(result_eigenvalue)


solve_neutralino_mass()
<<<0>>>