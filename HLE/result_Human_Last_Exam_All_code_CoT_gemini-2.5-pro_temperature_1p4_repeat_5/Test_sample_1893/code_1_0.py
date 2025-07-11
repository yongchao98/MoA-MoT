import sympy

def solve_neutralino_eigenvalue():
    """
    This function computes the specific eigenvalue of the neutralino mass matrix
    as requested by the user.
    """
    # Step 1 & 2: Define symbols and the matrix, then apply physical constraints.
    # We define the symbols for the parameters in the neutralino mass matrix.
    M, mu, Mz = sympy.symbols('M mu M_Z', real=True, positive=True)
    # The 'dynamic enhancement' condition requires M1=M2=M and cos(2*beta)=0, which implies sin(2*beta)=1.
    # The original matrix is:
    # M_N = sympy.Matrix([
    #     [M1*cW**2 + M2*sW**2, (M2 - M1)*sW*cW, 0, 0],
    #     [(M2 - M1)*sW*cW, M1*sW**2 + M2*cW**2, Mz, 0],
    #     [0, Mz, mu*sympy.sin(2*beta), -mu*sympy.cos(2*beta)],
    #     [0, 0, -mu*sympy.cos(2*beta), -mu*sympy.sin(2*beta)]
    # ])
    # Applying the constraints M1=M2=M simplifies M1*cW**2+M2*sW**2 to M*(cW**2+sW**2) = M.
    # The (1,2) and (2,1) elements become 0.
    # The constraint cos(2*beta)=0 makes the (3,4) and (4,3) elements 0.
    # The constraint sin(2*beta)=1 simplifies the other terms.
    M_N_simplified = sympy.Matrix([
        [M, 0, 0, 0],
        [0, M, Mz, 0],
        [0, Mz, mu, 0],
        [0, 0, 0, -mu]
    ])

    # Step 3: Determine the eigenvalues of the simplified matrix.
    # By inspection of the block-diagonal matrix, two eigenvalues are M and -mu.
    # These are proportional to the adjustable parameters M (which is M1=M2) and mu.
    # The other two eigenvalues are derived from the central 2x2 sub-matrix:
    # | M   Mz |
    # | Mz  mu |
    # Its characteristic equation is (M-lambda)(mu-lambda) - Mz**2 = 0.
    # The roots are lambda = (M+mu +/- sqrt((M-mu)**2 + 4*Mz**2)) / 2.
    # These two eigenvalues are not proportional to M or mu.

    other_eigs_expr = [
        (M + mu + sympy.sqrt((M - mu)**2 + 4*Mz**2)) / 2,
        (M + mu - sympy.sqrt((M - mu)**2 + 4*Mz**2)) / 2
    ]

    # Step 4: Compute the final value by considering a specific limit.
    # To 'compute' a value, we consider the physically motivated scenario where the
    # gamma-tilde and H_b-tilde states are light, meaning M << Mz and |mu| << Mz.
    # We find the value of these eigenvalues in the limit M -> 0 and mu -> 0.

    eig_limit_pos = sympy.limit(other_eigs_expr[0], M, 0)
    eig_limit_pos = sympy.limit(eig_limit_pos, mu, 0)

    eig_limit_neg = sympy.limit(other_eigs_expr[1], M, 0)
    eig_limit_neg = sympy.limit(eig_limit_neg, mu, 0)

    # The resulting eigenvalues are Mz and -Mz.
    # We will present the positive eigenvalue.
    final_eigenvalue_symbolic = eig_limit_pos
    
    # Use the known value for the Z boson mass in GeV.
    mz_value = 91.1876
    final_eigenvalue_numeric = final_eigenvalue_symbolic.subs(Mz, mz_value)

    # Print the result as requested.
    print("The eigenvalue (lambda) that is not proportional to M1, M2, or mu has been computed.")
    print("Under the specified conditions and in the limit where the adjustable mass parameters M and mu approach zero, the final equation for the eigenvalue is:")
    print(f"lambda = {final_eigenvalue_symbolic}")
    print("\nUsing the standard value for the Z boson mass (M_Z):")
    print(f"The number for M_Z is: {mz_value} GeV")
    print(f"The computed value of the eigenvalue is: {final_eigenvalue_numeric} GeV")

if __name__ == '__main__':
    solve_neutralino_eigenvalue()
<<<91.1876>>>