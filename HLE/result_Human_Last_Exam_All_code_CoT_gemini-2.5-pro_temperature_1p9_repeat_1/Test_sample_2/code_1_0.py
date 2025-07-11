import collections

def compute_spin_bordism_bg2_12():
    """
    Computes the reduced 12-th dimensional Spin bordism of BG2
    by analyzing the Atiyah-Hirzebruch Spectral Sequence.
    """
    # Step 1: Data for Spin bordism groups of a point Omega_q^{Spin}(pt).
    # Format: (rank of free part, list of torsion parts as strings).
    spin_bordism_point = {
        0: (1, []), 1: (0, ['Z2']), 2: (0, ['Z2']), 3: (0, []),
        4: (1, []), 5: (0, []), 6: (0, []), 7: (0, []),
        8: (1, []), 9: (0, ['Z2', 'Z2']), 10: (0, ['Z2']),
        11: (0, []), 12: (1, [])
    }

    # Step 2: Data for integral homology groups H_p(BG2; Z).
    # H*(BG2; Z) is torsion-free with H*(BG2; Q) = Q[x_4, x_12].
    # Ranks of H_p(BG2; Z) are determined by this ring structure.
    # rank H_p = number of ways to write p as 4a + 12b for integers a, b >= 0.
    homology_bg2_ranks = collections.defaultdict(int)
    homology_bg2_ranks[0] = 1 # from constant 1
    homology_bg2_ranks[4] = 1 # from x_4
    homology_bg2_ranks[8] = 1 # from x_4^2
    homology_bg2_ranks[12] = 2 # from x_4^3 and x_12
    
    total_degree = 12
    total_free_rank = 0
    total_torsion_parts = []
    
    print(f"Calculating the E^2 page terms E^2_{{p,q}} = H_p(BG2) \otimes Omega_q^{{Spin}} for p+q = {total_degree} (reduced case, p>0):")
    print("-" * 80)

    # We sum the ranks of the contributing terms E^2_{p,q}.
    # We will justify that all differentials are zero, so E^2 = E^infty.
    # And since all contributors are free, the group is their direct sum.
    
    components_ranks = []

    # Loop through p from 1 to total_degree for the reduced group
    for p in range(1, total_degree + 1):
        q = total_degree - p
        
        hp_rank = homology_bg2_ranks[p]
        if hp_rank == 0:
            continue
            
        omega_q_rank, omega_q_torsion = spin_bordism_point.get(q, (0, []))
        
        # Since H_*(BG2; Z) is free, E^2_{p,q} = H_p(BG2) tensor Omega_q
        # and H_p(BG2) = Z^hp_rank.
        term_rank = hp_rank * omega_q_rank
        # Torsion part of E^2 term: (Z/2)^hp_rank for Omega_q with Z/2, etc.
        term_torsion = omega_q_torsion * hp_rank
        
        if term_rank > 0 or term_torsion:
            components_ranks.append(term_rank)
            print(f"p={p}, q={q}:")
            print(f"  H_{p}(BG2; Z) = Z^{hp_rank}")
            omega_q_str = f"Z^{omega_q_rank}" if omega_q_rank > 0 else ""
            if omega_q_torsion:
                omega_q_str += " + " if omega_q_str else ""
                omega_q_str += " + ".join(omega_q_torsion)
            if not omega_q_str: omega_q_str = "0"
            print(f"  Omega_{q}^{{Spin}}(pt) = {omega_q_str}")

            e2_term_str = f"Z^{term_rank}" if term_rank > 0 else ""
            if term_torsion:
                e2_term_str += " + " if e2_term_str else ""
                e2_term_str += " + ".join(term_torsion)

            print(f"  => E^2_{{{p},{q}}} contributes a rank of {term_rank} to the final group.")
            print("-" * 80)
            
            total_free_rank += term_rank
            total_torsion_parts.extend(term_torsion)
            
    print("Summary of Calculation:")
    print("The Atiyah-Hirzebruch spectral sequence collapses at the E^2 page for total degree 12.")
    print("This is because the homology groups H_k(BG2; Z) are zero for k not in {0, 4, 8, 12, ...},")
    print("causing all potentially non-zero differentials d^r to have a zero source or target group.")
    print("\nThe contributing terms are all free abelian groups, so the final group is their direct sum.")

    equation = " + ".join(map(str, components_ranks))
    print(f"\nThe ranks of the components sum up: {equation} = {total_free_rank}")

    final_group_str = f"Z^{total_free_rank}" if total_free_rank > 0 else "0"
    if total_torsion_parts:
        final_group_str += " + " + " + ".join(total_torsion_parts)
    print(f"\nThus, the reduced 12-th Spin bordism of BG2 is: {final_group_str}")
    
compute_spin_bordism_bg2_12()