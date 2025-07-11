def solve_rank_of_kernel():
    """
    This function calculates the rank of the kernel of the Hurewicz homomorphism
    for the described topological space Y.
    """
    # The orders of the cyclic groups forming the fundamental group of Y.
    # pi_1(Y) = Z_n1 * Z_n2 * Z_n3
    n1 = 5
    n2 = 8
    n3 = 2
    print(f"The fundamental group of Y is the free product G = Z_{n1} * Z_{n2} * Z_{n3}.")

    # The Hurewicz map sends pi_1(Y) to H_1(Y) = pi_1(Y)^ab.
    # The kernel K is the commutator subgroup [G, G].
    # The index of K in G is the order of the abelianization H_1(Y).
    # H_1(Y) = Z_n1 + Z_n2 + Z_n3.
    index = n1 * n2 * n3
    print(f"The index of the kernel K in G is |H_1(Y)| = {n1} * {n2} * {n3} = {index}.")

    # The Euler characteristic of a group G = G1 * G2 * G3 is chi(G) = chi(G1) + chi(G2) + chi(G3) - 2.
    # For a finite group H, chi(H) = 1/|H|.
    # We calculate chi(G) using fractions to be exact.
    # chi_G = 1/n1 + 1/n2 + 1/n3 - 2
    # To avoid floating point arithmetic, we calculate chi_K directly using integers.
    # chi_K = index * chi_G = index * (1/n1 + 1/n2 + 1/n3 - 2)
    # chi_K = (n2*n3) + (n1*n3) + (n1*n2) - 2 * (n1*n2*n3)
    chi_k = (n2 * n3) + (n1 * n3) + (n1 * n2) - 2 * index
    print(f"The Euler characteristic of the kernel K is chi(K) = {index} * (1/{n1} + 1/{n2} + 1/{n3} - 2) = {chi_k}.")

    # The kernel K is a free group of rank r. The Euler characteristic of a free group
    # of rank r is chi(K) = 1 - r.
    # From this, we can find the rank r.
    rank = 1 - chi_k
    print("The rank r of K is determined by the equation chi(K) = 1 - r.")
    print(f"Substituting the value of chi(K), we get the equation:")
    print(f"{chi_k} = 1 - r")
    print(f"Solving for r yields:")
    print(f"r = 1 - ({chi_k}) = {rank}")

solve_rank_of_kernel()