import sympy

def compute_neutralino_eigenvalue():
    """
    This function symbolically computes the eigenvalues of the neutralino mass matrix
    under the specified physical conditions and identifies the one that is not
    proportional to the adjustable mass parameters M_1, M_2, or mu.
    """
    # Define the symbolic variables from the problem
    # M_Z is the mass of the Z boson, a known constant.
    M1, M2, mu, M_Z = sympy.symbols('M_1 M_2 mu M_Z', real=True, positive=True)
    beta, theta_W = sympy.symbols('beta theta_W', real=True)

    # For simplicity, we can directly implement the conditions for dynamic enhancement
    # which are M1=M2=M and beta=pi/4, leading to sin(2*beta)=1 and cos(2*beta)=0.
    M = sympy.Symbol('M', real=True, positive=True)

    # The simplified neutralino mass matrix under the enhancement conditions
    # As derived in the explanation, the matrix becomes block-diagonal.
    # The relevant block for the non-trivial eigenvalues mixes the ~Z and ~H_a states.
    sub_matrix = sympy.Matrix([
        [M, M_Z],
        [M_Z, mu]
    ])

    # Find the eigenvalues of this 2x2 sub-matrix
    eigenvalues = sub_matrix.eigenvals().keys()

    # The physical scenario requires the ~gamma and ~H_b states (with masses M and |mu|)
    # to be the lightest. This imposes the condition that M and mu must be small
    # compared to M_Z. We find the eigenvalues in the limit where M -> 0 and mu -> 0.
    limit_eigenvalues = []
    for ev in eigenvalues:
        # Sequentially take the limits
        ev_limit = sympy.limit(ev, M, 0)
        ev_limit = sympy.limit(ev_limit, mu, 0)
        limit_eigenvalues.append(ev_limit)

    # The eigenvalues in this limit are M_Z and -M_Z. They are not proportional
    # to M_1, M_2, or mu. We are asked for "the" eigenvalue.
    # We will output the positive one, which is the mass of the Z boson.
    final_eigenvalue = [val for val in limit_eigenvalues if val != 0][0]

    # The question requests to output each number in the final equation.
    # The final answer is symbolic, M_Z.
    print("The derived eigenvalue is:")
    # Using str() to print the symbol name
    print(str(final_eigenvalue))

compute_neutralino_eigenvalue()