import sympy

def solve_neutralino_mass():
    """
    Calculates the specific eigenvalue of the neutralino mass matrix
    under the condition of dynamic enhancement.
    """
    # Define the symbolic parameters from the neutralino mass matrix
    M1, M2, mu, beta, M_Z = sympy.symbols('M_1 M_2 mu beta M_Z')
    sw, cw = sympy.symbols('s_W c_W') # sin(theta_W), cos(theta_W)

    # The neutralino mass matrix in the specified basis
    # Basis: (-i*gamma_tilde, -i*Z_tilde, H_a_tilde, H_b_tilde)
    M_n = sympy.Matrix([
        [M1*cw**2 + M2*sw**2, (M2 - M1)*sw*cw, 0, 0],
        [(M2 - M1)*sw*cw, M1*sw**2 + M2*cw**2, M_Z, 0],
        [0, M_Z, mu*sympy.sin(2*beta), -mu*sympy.cos(2*beta)],
        [0, 0, -mu*sympy.cos(2*beta), -mu*sympy.sin(2*beta)]
    ])

    # The condition for dynamic enhancement is that the states gamma_tilde (1st basis vector)
    # and H_b_tilde (4th basis vector) do not mix with Z_tilde (2nd) and H_a_tilde (3rd).
    # This requires the off-diagonal blocks connecting the (1,4) and (2,3) subspaces to be zero.
    # M_n[0,1] = (M2 - M1)*sw*cw = 0  => M1 = M2 (since sw, cw != 0)
    # M_n[2,3] = -mu*cos(2*beta) = 0 => cos(2*beta) = 0 (assuming mu != 0)
    # For 0 <= beta <= pi/2, cos(2*beta) = 0 implies 2*beta = pi/2, so beta = pi/4.

    # Apply these conditions to the matrix
    M_n_simplified = M_n.subs({
        M2: M1,
        sympy.sin(2*beta): sympy.sin(sympy.pi/2), # sin(pi/2) = 1
        sympy.cos(2*beta): sympy.cos(sympy.pi/2)  # cos(pi/2) = 0
    })
    
    # The problem asks for the eigenvalue not proportional to M1, M2, or mu.
    # The eigenvalues of the simplified matrix still depend on M1 and mu.
    # To find a value independent of these adjustable parameters, we consider
    # the limit where they are zero.
    M_n_limit = M_n_simplified.subs({
        M1: 0,
        mu: 0
    })

    # Calculate the eigenvalues of this final matrix
    eigenvalues = M_n_limit.eigenvals() # Returns a dictionary {eigenvalue: multiplicity}

    # Print the results and the reasoning
    print("Step 1: The physical condition of dynamic enhancement simplifies the parameters to M_1 = M_2 and beta = pi/4.")
    print("Step 2: The simplified mass matrix is:")
    sympy.pprint(M_n_simplified)
    print("\nStep 3: To find an eigenvalue not proportional to the adjustable parameters M_1 or mu, we consider the limit where M_1 -> 0 and mu -> 0.")
    print("The matrix in this limit is:")
    sympy.pprint(M_n_limit)
    print("\nStep 4: The eigenvalues of this final matrix are:")
    print(list(eigenvalues.keys()))
    
    # Find the non-zero, positive eigenvalue from the list
    final_eigenvalue = None
    for val in eigenvalues.keys():
        if val > 0:
            final_eigenvalue = val
            break
            
    print(f"\nStep 5: The eigenvalue which is not proportional to M_1, M_2, or mu is the non-zero physical mass in this limit.")
    print(f"The final equation for the eigenvalue is lambda = {final_eigenvalue}.")
    # The prompt asks to output each number in the final equation. Here the equation is simple.
    # The numbers are implicitly 1.
    print(f"The value is 1 * {final_eigenvalue}.")


solve_neutralino_mass()