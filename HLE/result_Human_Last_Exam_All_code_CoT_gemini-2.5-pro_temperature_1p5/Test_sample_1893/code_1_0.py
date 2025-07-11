import sympy

def solve_neutralino_mass():
    """
    Computes the specific eigenvalue of the neutralino mass matrix based on the problem's conditions.
    """
    # Define symbols for the parameters in the matrix
    M1, M2, mu, beta, theta_W, M_Z = sympy.symbols('M1 M2 mu beta theta_W M_Z', real=True)
    
    # Define Weinberg angle functions and beta functions
    sW = sympy.sin(theta_W)
    cW = sympy.cos(theta_W)
    s2b = sympy.sin(2 * beta)
    c2b = sympy.cos(2 * beta)

    # Define the neutralino mass matrix as given in the problem
    M_N = sympy.Matrix([
        [M1 * cW**2 + M2 * sW**2, (M2 - M1) * sW * cW, 0, 0],
        [(M2 - M1) * sW * cW, M1 * sW**2 + M2 * cW**2, M_Z, 0],
        [0, M_Z, mu * s2b, -mu * c2b],
        [0, 0, -mu * c2b, -mu * s2b]
    ])

    # Apply the decoupling conditions for dynamic enhancement: M1=M2 and beta=pi/4
    # This makes the matrix block-diagonal.
    M_decoupled = M_N.subs({
        M2: M1,
        beta: sympy.pi / 4
    })

    # To find the eigenvalue not proportional to M1, M2, or mu, we consider the limit
    # where these parameters are zero.
    M_limit = M_decoupled.subs({
        M1: 0,
        mu: 0
    })

    # The matrix in this limit is:
    # [[0, 0,  0, 0],
    #  [0, 0, M_Z, 0],
    #  [0, M_Z,  0, 0],
    #  [0, 0,  0, 0]]
    
    # Calculate the eigenvalues of this final matrix.
    # The characteristic polynomial is L^2 * (L^2 - M_Z^2) = 0.
    # The eigenvalues are 0 (with multiplicity 2), M_Z, and -M_Z.
    eigenvalues = M_limit.eigenvals()

    # The physical value for the Z boson mass in GeV
    M_Z_value = 91.1876

    # Extract the non-zero eigenvalues and substitute the known value for M_Z
    result_eigenvalue = None
    for eig in eigenvalues.keys():
        if eig != 0:
            # We select the positive eigenvalue, which corresponds to the physical mass M_Z.
            if eig.is_positive:
                result_eigenvalue = eig.subs({M_Z: M_Z_value})

    # The problem statement requires printing the equation as well.
    # The characteristic equation for the non-trivial 2x2 block is det([[-L, M_Z],[M_Z, -L]])=0
    # which simplifies to L^2 - M_Z^2 = 0.
    L = sympy.symbols('L')
    final_eq_lhs = L**2 - M_Z**2
    
    print("The characteristic equation for the relevant sub-block is:")
    sympy.pprint(final_eq_lhs, use_unicode=False) 
    print("= 0")
    print(f"\nSubstituting M_Z = {M_Z_value}, the equation becomes:")
    print(f"L^2 - {M_Z_value**2} = 0")
    print(f"L^2 = {M_Z_value**2}")
    print("\nThe positive eigenvalue solution is:")
    print(result_eigenvalue)

solve_neutralino_mass()