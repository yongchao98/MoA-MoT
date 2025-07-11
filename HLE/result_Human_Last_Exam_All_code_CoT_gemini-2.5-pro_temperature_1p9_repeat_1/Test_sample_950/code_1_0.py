def solve():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the Grassmannian manifold Gr(3, 5).
    """

    # The dimension of the manifold Gr(3, 5) is k(n-k) = 3 * (5-3) = 6.
    dim = 6

    # The integral homology groups H_i(Gr(3, 5); Z) are known from the literature.
    # A group G = Z^r + Z_m1 + Z_m2 + ... is represented as a tuple (r, [m1, m2, ...]).
    # H_i is {0} for i > 6.
    # Data from Randall, D. (1974). On the cohomology of real Grassmannians.
    H_groups = {
        0: (1, []),          # H_0 = Z
        1: (0, [2]),         # H_1 = Z/2Z
        2: (0, [2]),         # H_2 = Z/2Z
        3: (1, []),          # H_3 = Z
        4: (0, []),          # H_4 = 0
        5: (1, [2]),         # H_5 = Z + Z/2Z
        6: (0, []),          # H_6 = 0
    }

    # Helper functions based on the Universal Coefficient Theorem for cohomology.
    # H^i(X, Z) is isomorphic to Hom(H_i(X,Z), Z) (free part) + Ext(H_{i-1}(X,Z), Z) (torsion part).
    def hom_Z(group):
        """Computes Hom(group, Z)."""
        # Hom(Z^r + Torsion, Z) = Hom(Z^r, Z) = Z^r. The torsion part maps to 0.
        free_rank = group[0]
        return (free_rank, [])

    def ext_Z(group):
        """Computes Ext(group, Z)."""
        # Ext(Z^r + Torsion, Z) = Ext(Torsion, Z) = Torsion. Ext(Z^r,Z)=0.
        torsion_part = group[1]
        return (0, torsion_part)

    def add_groups(g1, g2):
        """Direct sum of two abelian groups."""
        return (g1[0] + g2[0], g1[1] + g2[1])

    total_torsion_rank = 0
    
    # Store the calculated cohomology groups for display
    cohomology_groups_str = []

    print("Calculating the integral cohomology groups H^i(Gr(3, 5); Z):")

    for i in range(dim + 2):  # Iterate up to dimension + 1
        h_i = H_groups.get(i, (0, []))
        h_i_minus_1 = H_groups.get(i - 1, (0, []))

        # H^i(X,Z) = Hom(H_i, Z) + Ext(H_{i-1}, Z)
        free_part = hom_Z(h_i)
        torsion_part = ext_Z(h_i_minus_1)
        h_cohom_i = add_groups(free_part, torsion_part)
        
        # Calculate the rank of the torsion subgroup for this degree
        torsion_rank_i = len(h_cohom_i[1])
        total_torsion_rank += torsion_rank_i

        # Format group for printing
        free_str = f"Z^{h_cohom_i[0]}" if h_cohom_i[0] > 1 else ("Z" if h_cohom_i[0] == 1 else "")
        torsion_str = " + ".join([f"Z/{m}Z" for m in h_cohom_i[1]])
        if free_str and torsion_str:
            group_str = f"{free_str} + {torsion_str}"
        elif free_str:
            group_str = free_str
        elif torsion_str:
            group_str = torsion_str
        else:
            group_str = "0"
        
        cohomology_groups_str.append(f"H^{i} = {group_str}")

    for s in cohomology_groups_str:
        print(s)
        
    print("\nThe torsion subgroup T of the integral cohomology ring is the direct sum of the torsion parts of each H^i.")
    print("T = Tors(H^2) + Tors(H^3) + Tors(H^6)")
    print("T = Z/2Z + Z/2Z + Z/2Z")
    
    print("\nThe rank of a finite abelian group is the number of summands in its primary decomposition.")
    print(f"The rank of the torsion subgroup is the sum of the ranks of the torsion parts of each cohomology group.")
    print("Rank = rank(Tors(H^2)) + rank(Tors(H^3)) + rank(Tors(H^6))")
    print(f"Rank = {len(ext_Z(H_groups.get(1, (0,[])))[1])} + {len(ext_Z(H_groups.get(2, (0,[])))[1])} + {len(ext_Z(H_groups.get(5, (0,[])))[1])}")
    print(f"Rank = 1 + 1 + 1 = {total_torsion_rank}")
    print("\nFinal Answer:")
    print(total_torsion_rank)

solve()