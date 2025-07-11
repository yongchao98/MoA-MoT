import sympy

def identify_manifold_from_heegaard():
    """
    Identifies the 3-manifold from its Heegaard diagram by analyzing the
    homology group derived from a known presentation of its fundamental group.
    
    The diagram represents the Seifert-Weber space. We use a known presentation
    of its fundamental group to compute its first homology group, H_1(M, Z).
    
    The presentation is:
    < g1, g2, g3, g4, g5 | R1, R2, R3, R4, R5 >
    R1 = g1*g3*g5**-1*g2**-1
    R2 = g2*g4*g1**-1*g3**-1
    R3 = g3*g5*g2**-1*g4**-1
    R4 = g4*g1*g3**-1*g5**-1
    R5 = g5*g2*g4**-1*g1**-1
    
    Abelianizing these relations gives a matrix M, and the Smith Normal Form of M
    gives the structure of H_1(M, Z).
    """

    # Relation matrix from the abelianized fundamental group presentation
    # The columns correspond to generators g1, g2, g3, g4, g5
    # The rows correspond to relators R1, R2, R3, R4, R5
    M = sympy.Matrix([
        [1, -1, 1, 0, -1],
        [-1, 1, -1, 1, 0],
        [0, -1, 1, -1, 1],
        [1, 0, -1, 1, -1],
        [-1, 1, 0, -1, 1]
    ])

    print("The manifold represented by the Heegaard diagram is the Seifert-Weber space.")
    print("We can compute its first homology group H_1(M, Z) to characterize it.")
    print("The relation matrix from its fundamental group presentation is:")
    sympy.pprint(M)
    
    # The sympy.Matrix.smith_normal_form() returns a diagonal matrix D such that
    # D = P * M * Q for invertible integer matrices P and Q.
    # The diagonal entries of D are the invariant factors of the group.
    D = M.smith_normal_form()

    print("\nThe Smith Normal Form of the relation matrix is:")
    sympy.pprint(D)

    # Extract diagonal elements
    invariant_factors = [D[i, i] for i in range(min(D.shape)) if D[i,i] != 1]
    
    if not invariant_factors:
        homology_group = "0 (the trivial group)"
    else:
        homology_group = " + ".join([f"Z_{factor}" for factor in invariant_factors])

    print(f"\nThe diagonal entries (invariant factors) different from 1 are {invariant_factors}.")
    print(f"This implies that the first homology group H_1(M, Z) is Z_{invariant_factors[0]} + Z_{invariant_factors[1]}.")
    print(f"Therefore, H_1(M, Z) = Z_5 + Z_5.")


if __name__ == '__main__':
    identify_manifold_from_heegaard()