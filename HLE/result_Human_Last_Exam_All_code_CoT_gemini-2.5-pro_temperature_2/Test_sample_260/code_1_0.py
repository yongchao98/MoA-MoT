from fractions import Fraction

def solve_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # The orders of the fundamental groups of X1, X2, and X3
    n1 = 5  # |pi_1(X1)| = |Z_5|
    n2 = 8  # |pi_1(X2)| = |Z_8|
    n3 = 2  # |pi_1(X3)| = |Z_2|
    
    # The number of spaces being connected summed
    k = 3
    
    orders = [n1, n2, n3]
    
    # G = Z_n1 * Z_n2 * Z_n3
    # K = Ker(h_*) is the commutator subgroup [G,G].
    # [G:K] = |H_1(Y)| = |Z_n1 + Z_n2 + Z_n3| = n1 * n2 * n3
    index_K = n1 * n2 * n3
    
    # The rational Euler characteristic of G is chi(G) = sum(1/|G_i|) - (k-1)
    chi_G_sum_part = sum(Fraction(1, n) for n in orders)
    chi_G = chi_G_sum_part - (k - 1)
    
    # For a subgroup K of G of finite index, chi(K) = [G:K] * chi(G).
    chi_K = index_K * chi_G
    
    # For a free group F_r of rank r, chi(F_r) = 1 - r.
    # Since K is a free group, chi(K) = 1 - rank(K).
    # So, rank(K) = 1 - chi(K).
    rank_K = 1 - chi_K
    
    print("Step 1: The fundamental groups are G1=Z_5, G2=Z_8, G3=Z_2.")
    print(f"The orders of these groups are |G1|={n1}, |G2|={n2}, |G3|={n3}.")
    
    print("\nStep 2: The fundamental group of the connected sum Y is G = G1 * G2 * G3.")
    print(f"The number of groups in the free product is k={k}.")
    
    print("\nStep 3: The kernel K of the Hurewicz map is the commutator subgroup of G.")
    print(f"The index of K in G is [G:K] = |H_1(Y)| = |G1|*|G2|*|G3| = {n1}*{n2}*{n3} = {index_K}.")

    print("\nStep 4: Calculate the rank of K using Euler characteristics.")
    print(f"The rational Euler characteristic of G is chi(G) = (1/{n1} + 1/{n2} + 1/{n3}) - ({k}-1) = {chi_G_sum_part} - {k-1} = {chi_G}.")
    print(f"The Euler characteristic of K is chi(K) = [G:K] * chi(G) = {index_K} * ({chi_G}) = {chi_K}.")
    print(f"The rank r of K is given by chi(K) = 1 - r.")
    print(f"Therefore, r = 1 - chi(K) = 1 - ({chi_K}) = {rank_K}.")
    
    print(f"\nThe final calculated rank is: {int(rank_K)}")

solve_rank()