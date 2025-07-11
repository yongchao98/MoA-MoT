def get_cohomology_torsion_rank(i, homology_groups):
    """
    Calculates the rank of the torsion subgroup of the i-th cohomology group
    using the Universal Coefficient Theorem.
    Tors(H^i) is isomorphic to Tors(H_{i-1}).
    All torsion in the homology of real Grassmannians is 2-torsion.
    The rank of (Z_2)^k is k.
    """
    # Get the (i-1)-th homology group. H[-1] is defined as 0.
    h_im1 = homology_groups.get(i - 1, {'free_rank': 0, 'torsion_rank': 0})
    
    # The torsion rank of H^i is the torsion rank of H_{i-1}
    return h_im1['torsion_rank']

def solve():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the Grassmannian Gr(3, 5).
    """
    # Homology groups H_i(Gr(3,5); Z) based on published tables.
    # 'torsion_rank' for (Z_2)^k is k.
    homology_groups = {
        0: {'free_rank': 1, 'torsion_rank': 0},      # H_0 = Z
        1: {'free_rank': 0, 'torsion_rank': 1},      # H_1 = Z_2
        2: {'free_rank': 1, 'torsion_rank': 0},      # H_2 = Z
        3: {'free_rank': 0, 'torsion_rank': 1},      # H_3 = Z_2
        4: {'free_rank': 1, 'torsion_rank': 1},      # H_4 = Z + Z_2
        5: {'free_rank': 0, 'torsion_rank': 2},      # H_5 = Z_2 + Z_2
        6: {'free_rank': 0, 'torsion_rank': 0},      # H_6 = 0
    }
    
    dim = 6 # Dimension of Gr(3,5)
    total_torsion_rank = 0
    torsion_ranks_per_degree = []

    print("Calculating the torsion rank for each cohomology group H^i(Gr(3,5); Z):")
    # We need to consider H^i up to dimension + 1 due to the UCT shift.
    for i in range(dim + 2):
        rank = get_cohomology_torsion_rank(i, homology_groups)
        print(f"  - Rank of Tors(H^{i}) = {rank}")
        if rank > 0:
            torsion_ranks_per_degree.append(str(rank))
        total_torsion_rank += rank
        
    equation = " + ".join(torsion_ranks_per_degree)
    
    print("\nThe rank of the torsion subgroup of the integral cohomology ring is the sum of these ranks.")
    print(f"Final calculation: {equation} = {total_torsion_rank}")
    print("\nThus, the total rank of the torsion subgroup is:")
    print(total_torsion_rank)

solve()