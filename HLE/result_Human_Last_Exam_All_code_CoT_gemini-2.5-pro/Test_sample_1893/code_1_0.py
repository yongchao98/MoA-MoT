import sympy

def solve_neutralino_mass():
    """
    Computes the special eigenvalue of the neutralino mass matrix under the
    condition of dynamic enhancement of radiative decay.

    The condition implies that some parameters must take specific values, leading
    to an eigenvalue that is independent of the other adjustable parameters.
    """
    # Define symbolic variables
    M1, M2, mu, beta, MZ = sympy.symbols('M1 M2 mu beta MZ')
    cW, sW = sympy.symbols('cW sW') # cos(theta_W), sin(theta_W)

    # The neutralino mass matrix from the problem description
    M_N = sympy.Matrix([
        [M1*cW**2 + M2*sW**2, (M2 - M1)*sW*cW, 0, 0],
        [(M2 - M1)*sW*cW, M1*sW**2 + M2*cW**2, MZ, 0],
        [0, MZ, mu*sympy.sin(2*beta), -mu*sympy.cos(2*beta)],
        [0, 0, -mu*sympy.cos(2*beta), -mu*sympy.sin(2*beta)]
    ])

    # The "dynamic enhancement" requires the lightest neutralinos to be pure
    # states. This implies the decoupling of certain states, which mathematically
    # means the mixing terms are zero. One way this is achieved is if mu = 0.
    # Let's set mu = 0 and see its effect on the eigenvalues.
    
    print("The condition of 'dynamic enhancement' can be interpreted as requiring the Higgsino mass parameter mu to be zero.")
    print("Let's analyze the matrix for mu = 0:")
    
    # Substitute mu = 0 into the matrix
    M_N_mu0 = M_N.subs(mu, 0)
    
    print("\nThe Neutralino Mass Matrix with mu = 0:")
    sympy.pprint(M_N_mu0)
    
    # A matrix with a zero row or column must have a determinant of 0.
    # This implies that one of its eigenvalues must be 0.
    
    # Let's verify this by calculating the determinant.
    det_mu0 = M_N_mu0.det()
    
    print("\nThe determinant of this matrix is:")
    sympy.pprint(det_mu0)
    
    print("\nSince the determinant is 0, one of the eigenvalues must be 0.")
    print("This eigenvalue is independent of the other parameters M1, M2, and MZ.")
    
    # The final equation is lambda = 0. We output the single number in this equation.
    final_eigenvalue = 0
    print(f"\nThe eigenvalue not proportional to M1, M2, or mu is:")
    print(final_eigenvalue)


solve_neutralino_mass()