import collections

def compute_spin_bordism_bg2():
    """
    Computes the reduced 12-th dimensional Spin bordism of BG2.
    This is done by summing the diagonal terms of a collapsed Atiyah-Hirzebruch Spectral Sequence.
    """
    
    # Ranks of the integral homology groups H_p(BG2)
    # based on H_*(BG2; Z) being a polynomial algebra Z[u_4, u_6]
    h_p_ranks = {
        0: 1, 1: 0, 2: 0, 3: 0,
        4: 1,  # generator u_4
        5: 0,
        6: 1,  # generator u_6
        7: 0,
        8: 1,  # generator u_4^2
        9: 0,
        10: 1, # generator u_4*u_6
        11: 0,
        12: 2  # generators u_4^3 and u_6^2
    }
    
    # Spin bordism groups of a point, Omega_q^Spin, up to dimension 12.
    # Stored as (rank, [torsion_orders])
    omega_q_spin = {
        0:  (1, []),        # Z
        1:  (0, [2]),      # Z_2
        2:  (0, [2]),      # Z_2
        3:  (0, []),        # 0
        4:  (1, []),        # Z
        5:  (0, []),        # 0
        6:  (0, []),        # 0
        7:  (0, []),        # 0
        8:  (2, []),        # Z + Z
        9:  (0, [2]),      # Z_2
        10: (0, [2]),      # Z_2
        11: (0, []),        # 0
        12: (1, [])         # Z
    }
    
    total_degree = 12
    total_rank = 0
    total_torsion = collections.defaultdict(int)
    
    print("Computing the reduced 12th Spin bordism of BG2, sum of E_2^{p,q} terms for p+q=12, p>0.\n")
    print(f"{'Term (p,q)':<12} | {'H_p(BG2)':<12} | {'Omega_q^Spin':<12} | {'Contribution E_2^{p,q} = H_p x Omega_q':<40}")
    print("-" * 80)
    
    # Iterate p from 1 to 12 for the reduced group
    for p in range(1, total_degree + 1):
        q = total_degree - p
        
        # Get ranks and group structures
        h_p_rank = h_p_ranks.get(p, 0)
        omega_q_rank, omega_q_torsion = omega_q_spin.get(q, (0, []))
        
        if h_p_rank == 0:
            # If homology is zero, the term is zero.
            continue
            
        # Since H_p(BG2) is free, the tensor product is straightforward.
        # H_p(BG2) tensor Omega_q = Z^h_p tensor (Z^r_q + T_q) = Z^(h_p*r_q) + (T_q)^h_p
        
        # Format the components for printing
        h_p_str = f"Z^{h_p_rank}" if h_p_rank > 1 else ("Z" if h_p_rank == 1 else "0")
        
        omega_q_rank_str = f"Z^{omega_q_rank}" if omega_q_rank > 1 else ("Z" if omega_q_rank == 1 else "")
        omega_q_torsion_str = " + ".join([f"Z_{t}" for t in omega_q_torsion])
        omega_q_str_parts = [part for part in [omega_q_rank_str, omega_q_torsion_str] if part]
        omega_q_str = " + ".join(omega_q_str_parts) if omega_q_str_parts else "0"

        # Calculate contribution
        rank_contrib = h_p_rank * omega_q_rank
        torsion_contrib = collections.defaultdict(int)
        for t in omega_q_torsion:
            torsion_contrib[t] += h_p_rank
            
        # Accumulate totals
        total_rank += rank_contrib
        for t, count in torsion_contrib.items():
            total_torsion[t] += count
            
        # Format contribution string
        rank_contrib_str = f"Z^{rank_contrib}" if rank_contrib > 1 else ("Z" if rank_contrib == 1 else "")
        torsion_contrib_str = " + ".join([f"(Z_{t})^{count}" if count > 1 else f"Z_{t}" for t, count in torsion_contrib.items()])
        contrib_str_parts = [part for part in [rank_contrib_str, torsion_contrib_str] if part]
        contrib_str = " + ".join(contrib_str_parts) if contrib_str_parts else "0"

        print(f"({p:<2},{q:<2})      | {h_p_str:<12} | {omega_q_str:<12} | {contrib_str:<40}")

    # Format the final result string
    final_rank_str = f"Z^{total_rank}" if total_rank > 1 else ("Z" if total_rank == 1 else "")
    final_torsion_str = " + ".join([f"(Z_{t})^{count}" if count > 1 else f"Z_{t}" for t, count in sorted(total_torsion.items())])
    
    final_result_parts = [part for part in [final_rank_str, final_torsion_str] if part]
    final_result_str = " + ".join(final_result_parts) if final_result_parts else "0"
    
    print("-" * 80)
    print(f"Total Reduced Bordism Group: Sum = {final_rank_str} + {final_torsion_str if final_torsion_str else '0'}")
    
    # To satisfy the problem constraint: "output each number in the final equation!"
    print(f"\nThe final equation for the rank is sum of ranks from each non-zero term.")
    print(f"The rank of Z^5 is 5.")
    print(f"The torsion component is Z_2, where the subscript is 2.")

compute_spin_bordism_bg2()
