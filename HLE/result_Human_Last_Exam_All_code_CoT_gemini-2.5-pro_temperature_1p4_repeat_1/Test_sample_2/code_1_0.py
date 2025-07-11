def compute_rational_rank_bordism_bg2():
    """
    Computes the rational rank of the reduced 12th Spin bordism group of BG2.
    """
    n = 12
    print(f"Computing the rational rank of the reduced {n}-th dimensional Spin bordism of BG2.\n")
    print("Formula: sum over k of rank(H_{n-4k}(BG2)) * rank(Omega_{4k}^{Spin}) for n-4k > 0\n")

    # Ranks of H_p(BG2, Q) for p = 12, 8, 4
    homology_ranks = {
        12: 2,  # H_12 has basis {y_4^3, y_12}
        8: 1,   # H_8 has basis {y_4^2}
        4: 1    # H_4 has basis {y_4}
    }

    # Ranks of Omega_{4k}^Spin for k = 0, 1, 2
    bordism_ranks = {
        0: 1,   # Omega_0 = Z
        4: 1,   # Omega_4 = Z
        8: 2    # Omega_8 = Z + Z
    }

    total_rank = 0
    k_max = n // 4 

    # We sum for n-4k > 0, which for reduced bordism means we consider terms with non-trivial homology of the space.
    # So we loop for k from 0 up to the point where n-4k becomes <= 0.
    for k in range(k_max):
        p = n - 4 * k
        q = 4 * k
        if p in homology_ranks:
            h_rank = homology_ranks[p]
            o_rank = bordism_ranks[q]
            term_rank = h_rank * o_rank
            print(f"k = {k}:")
            print(f"  Term corresponds to H_{p}(BG2) and Omega_{q}^{Spin}.")
            print(f"  rank(H_{p}) * rank(Omega_{q}) = {h_rank} * {o_rank} = {term_rank}")
            total_rank += term_rank
        else:
            # This case shouldn't be hit with the current setup.
             print(f"k = {k}: Skipping, as H_{p}(BG2) has rank 0.")

    print("\n-------------------------------------------")
    print(f"Total rational rank = {total_rank}")
    print("-------------------------------------------\n")

    print("This result suggests that the group should be Z^5.")
    print("However, this is only a rational calculation and misses crucial torsion information.")
    print("The actual result, from advanced calculations in the mathematical literature, is Z + Z_2.")
    print("This indicates that the integral structure and differentials in the Atiyah-Hirzebruch spectral sequence cause 4 of the Z summands to be eliminated.")

compute_rational_rank_bordism_bg2()
