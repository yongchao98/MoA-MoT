import sympy

def solve_neutralino_eigenvalue():
    """
    This function computes the eigenvalues of the neutralino mass matrix
    under the condition of dynamic enhancement and identifies the one
    that is not proportional to M1, M2, or mu.
    """
    # Define symbolic variables for the matrix parameters
    M1, M2, mu, M_Z = sympy.symbols('M1 M2 mu M_Z', real=True, positive=True)
    theta_W, beta = sympy.symbols('theta_W beta', real=True)

    # Define sine and cosine of the Weinberg angle for brevity
    s_W = sympy.sin(theta_W)
    c_W = sympy.cos(theta_W)

    # Define the neutralino mass matrix
    M_N = sympy.Matrix([
        [M1*c_W**2 + M2*s_W**2, (M2 - M1)*s_W*c_W, 0, 0],
        [(M2 - M1)*s_W*c_W, M1*s_W**2 + M2*c_W**2, M_Z, 0],
        [0, M_Z, mu*sympy.sin(2*beta), -mu*sympy.cos(2*beta)],
        [0, 0, -mu*sympy.cos(2*beta), -mu*sympy.sin(2*beta)]
    ])

    # The condition for dynamic enhancement is that gamma_tilde and H_b_tilde are pure states.
    # This requires the off-diagonal elements connecting them to other states to be zero.
    # M_N[0, 1] = 0  => (M2 - M1)*s_W*c_W = 0  => M1 = M2 (since s_W, c_W != 0)
    # M_N[2, 3] = 0  => -mu*cos(2*beta) = 0    => mu = 0 or cos(2*beta) = 0
    #
    # The problem asks for an eigenvalue NOT proportional to M1, M2, or mu.
    # The case mu=0 provides such an eigenvalue, while the case cos(2*beta)=0 does not.
    # Therefore, we apply the conditions M2 = M1 and mu = 0.
    
    print("Applying conditions for dynamic enhancement: M2 = M1 and mu = 0")
    M_N_simplified = M_N.subs([(M2, M1), (mu, 0)])

    print("\nThe simplified neutralino mass matrix is:")
    sympy.pprint(M_N_simplified)

    # Calculate the eigenvalues of the simplified matrix.
    # The eigenvals() method returns a dictionary of {eigenvalue: multiplicity}.
    eigenvalues_dict = M_N_simplified.eigenvals()
    eigenvalues = list(eigenvalues_dict.keys())

    print("\nThe eigenvalues of this matrix are:")
    for eig in eigenvalues:
        sympy.pprint(eig)
        print("-" * 20)

    # Find the eigenvalue that is independent of M1 and mu.
    result_eigenvalue = None
    for eig in eigenvalues:
        # Check if the eigenvalue expression contains any of the variable parameters.
        if not eig.free_symbols.intersection({M1, M2, mu}):
            result_eigenvalue = eig
            break
            
    print("\nThe eigenvalue which is not proportional to M1, M2, or mu is:")
    print(result_eigenvalue)

if __name__ == '__main__':
    solve_neutralino_eigenvalue()