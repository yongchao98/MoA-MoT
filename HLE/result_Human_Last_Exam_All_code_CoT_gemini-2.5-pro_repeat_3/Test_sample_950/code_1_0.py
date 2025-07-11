import collections

def get_group_name(rank, torsion_coeffs):
    """Generates a string representation of a finitely generated abelian group."""
    parts = []
    if rank > 0:
        if rank == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{rank}")
    
    # Use a counter for prettier output like Z_2^2 instead of Z_2 + Z_2
    counts = collections.Counter(torsion_coeffs)
    for coeff, num in sorted(counts.items()):
        if num == 1:
            parts.append(f"Z_{coeff}")
        else:
            parts.append(f"(Z_{coeff})^{num}")
            
    if not parts:
        return "0"
    return " + ".join(parts)

def main():
    """
    Calculates the rank of the torsion subgroup of H*(G(3,5); Z).
    """
    print("The space of 3-subspaces in R^5 is the Grassmannian G(3,5).")
    print("This is diffeomorphic to G(2,5), a 6-dimensional manifold.")
    print("We will compute the rank of the torsion subgroup of its integral cohomology ring H*(G(2,5); Z).\n")

    # Step 1: Define the integral homology groups H_i(G(2,5); Z)
    # Each group is represented as a tuple (rank_of_free_part, list_of_torsion_coefficients)
    homology_groups = {
        -1: (0, []), # H_{-1} is trivial
        0: (1, []),   # H_0 = Z
        1: (0, [2]),  # H_1 = Z_2
        2: (0, [2]),  # H_2 = Z_2
        3: (0, []),   # H_3 = 0
        4: (1, []),   # H_4 = Z
        5: (0, []),   # H_5 = 0
        6: (1, []),   # H_6 = Z
    }

    print("Step 1: Integral Homology Groups of G(2,5)")
    for i in range(7):
        rank, torsion = homology_groups[i]
        print(f"H_{i}: {get_group_name(rank, torsion)}")
    print("-" * 30)

    # Step 2: Compute cohomology groups H^k using UCT
    cohomology_groups = {}
    total_torsion_coeffs = []
    
    print("Step 2: Calculating Integral Cohomology Groups H^k using UCT")
    print("H^k = Free(H_k) + Torsion(H_{k-1})\n")

    for k in range(7):
        # H^k = Free(H_k) + Torsion(H_{k-1})
        free_rank = homology_groups[k][0]
        torsion_part = homology_groups[k-1][1]
        
        cohomology_groups[k] = (free_rank, torsion_part)
        total_torsion_coeffs.extend(torsion_part)
        
        print(f"H^{k}: Free part from H_{k} (rank={homology_groups[k][0]}), Torsion part from H_{k-1} ({get_group_name(0, homology_groups[k-1][1])})")
        print(f"-> H^{k} = {get_group_name(free_rank, torsion_part)}\n")
    print("-" * 30)
    
    # Step 3: Find the total torsion subgroup and its rank
    print("Step 3: Analyzing the Torsion Subgroup")
    torsion_subgroup_str = get_group_name(0, total_torsion_coeffs)
    if not total_torsion_coeffs:
        torsion_subgroup_str = "0"
        
    print(f"The torsion subgroup of the integral cohomology ring is the direct sum of the torsion parts of each H^k.")
    
    torsion_components = []
    for k in range(7):
        torsion_part = cohomology_groups[k][1]
        if torsion_part:
            torsion_components.append(f"Torsion(H^{k}) = {get_group_name(0, torsion_part)}")
            
    if torsion_components:
        print(" + ".join(torsion_components), end=" = ")
    
    print(f"{torsion_subgroup_str}")

    # The rank of a finite abelian group is the minimum number of generators.
    # For a direct sum of cyclic groups Z_m1 + Z_m2 + ..., this is the number of summands.
    rank_of_torsion = len(total_torsion_coeffs)
    
    print(f"\nThe rank of this group (its minimal number of generators) is the number of its cyclic components.")
    print(f"The number of cyclic components in {torsion_subgroup_str} is {rank_of_torsion}.")
    print("-" * 30)
    print(f"Final Answer: The rank of the torsion subgroup is {rank_of_torsion}.")

if __name__ == '__main__':
    main()