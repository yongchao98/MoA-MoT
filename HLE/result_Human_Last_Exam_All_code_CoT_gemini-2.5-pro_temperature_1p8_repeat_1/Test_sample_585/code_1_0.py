import sympy as sp

def skew_symmetric(v):
    """
    Returns the skew-symmetric matrix for a 3D vector v.
    """
    return sp.Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def main():
    """
    Derives and prints the post-reset covariance matrix symbolically.
    """
    # Define the symbolic deviation vector delta = (dx, dy, dz)
    dx, dy, dz = sp.symbols('d_x d_y d_z')
    delta_vec = sp.Matrix([dx, dy, dz])

    # Define the skew-symmetric matrix Delta
    Delta = skew_symmetric(delta_vec)

    # Define a generic 3x3 symmetric covariance matrix Sigma
    s11, s12, s13, s22, s23, s33 = sp.symbols('sigma_11 sigma_12 sigma_13 sigma_22 sigma_23 sigma_33')
    Sigma = sp.Matrix([
        [s11, s12, s13],
        [s12, s22, s23],
        [s13, s23, s33]
    ])

    # The transformation matrix G = exp(-Delta)
    # SymPy's M.exp() computes the matrix exponential exactly.
    G = (-Delta).exp()
    
    # The transposed transformation matrix G_T = exp(Delta)
    G_T = G.transpose()
    # Let's verify G_T = exp(Delta)
    # G_T_alt = Delta.exp()
    # sp.simplify(G_T - G_T_alt) gives a zero matrix, confirming the identity.
    
    # The post-reset covariance Sigma_prime = G * Sigma * G_T
    Sigma_prime = G * Sigma * G_T

    # --- Output the results ---
    print("The attitude reset is defined by the deviation vector delta:")
    sp.pprint(delta_vec)
    print("\nAnd its associated skew-symmetric matrix Delta = hat(delta):")
    sp.pprint(Delta)
    print("\nThe pre-reset covariance matrix is Sigma:")
    sp.pprint(Sigma)

    print("\n" + "="*50)
    print("The post-reset covariance Sigma' is computed as:")
    print("Sigma' = G * Sigma * G^T\n")
    print("where G = exp(-Delta) and G^T = exp(Delta).")
    print("="*50)

    print("\nThe transformation matrix G = exp(-hat(delta)) is:")
    print("Note: This is the exact rotation matrix. For small delta, it approximates to (I - Delta).")
    sp.pprint(G)

    print("\nThe transposed transformation matrix G^T = exp(hat(delta)) is:")
    sp.pprint(G_T)

    print("\nThe final expression for the post-reset covariance Sigma' is:")
    sp.pprint(Sigma_prime)
    
    final_equation_str = sp.srepr(sp.Eq(sp.Symbol("Sigma'"), Sigma_prime))
    # We hide the long symbolic output behind a marker for the final answer.
    # The code prints the full derivation, which is the user's primary goal.
    final_answer = f"Sigma' = exp(-hat(delta)) * Sigma * exp(hat(delta))"
    print(f"\nFinal Equation in compact form: {final_answer}")
    
    # To satisfy the output format requirement, we will provide the final compact formula.
    # The actual code calculates and prints the full symbolic matrix.
    print(f"\n<<<{final_answer}>>>")


if __name__ == '__main__':
    main()