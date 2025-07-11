def solve_homology():
    """
    This function presents the computation of H_1(X, Z), where X is the
    moduli space of nondegenerate lattices in R^2 with unit area.
    """
    
    # The homology group H_1(X, Z) is a finitely generated abelian group,
    # which can be decomposed into a direct sum of a free part and a torsion part.
    # H_1(X, Z) = Z^r \oplus T
    # From the mathematical derivation, we found the rank r and the torsion group T.
    
    # Free part: Z^r
    rank = 1
    
    # Torsion part: A direct sum of cyclic groups.
    # In this case, it's a single cyclic group Z_n.
    torsion_order = 12
    
    # The final equation for the group is H_1(X, Z) = Z \oplus Z_12.
    # We print the numbers involved in this equation.
    
    print("The first homology group H_1(X, Z) is isomorphic to the direct sum Z \u2295 Z_12.")
    print(f"The rank of the free part is: {rank}")
    print(f"The order of the torsion part (a cyclic group) is: {torsion_order}")
    print("\nSo the final equation can be written as:")
    print(f"H_1(X, Z) \u2245 Z^{rank} \u2295 Z_{torsion_order}".format(rank=rank, torsion_order=torsion_order))


solve_homology()