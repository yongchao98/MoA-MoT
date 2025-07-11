def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Step 1: Define the orders of the finite groups from the spaces X_i
    n1 = 5  # Order of pi_1(X_1) = Z/5Z
    n2 = 8  # Order of pi_1(X_2) = Z/8Z
    n3 = 2  # Order of pi_1(X_3) = Z/2Z
    orders = [n1, n2, n3]
    num_groups = len(orders)

    print(f"The fundamental group of Y is G = Z/{n1}Z * Z/{n2}Z * Z/{n3}Z.")

    # Step 2: Calculate the index of the kernel K in G.
    # The kernel K is the commutator subgroup [G,G].
    # The index [G:K] is the order of the abelianization G/[G,G] = H_1(Y).
    index = 1
    for n in orders:
        index *= n
    
    print(f"The kernel K is the commutator subgroup [G,G].")
    print(f"The index of K in G is the order of H_1(Y) = Z/{n1}Z + Z/{n2}Z + Z/{n3}Z.")
    print(f"Index [G:K] = {n1} * {n2} * {n3} = {index}.")

    # Step 3: Calculate the rational Euler characteristic of G.
    # For a finite group H, chi(H) = 1/|H|.
    # For a free product G = G1*G2*...*Gn, chi(G) = sum(chi(Gi)) - (n-1).
    chi_G_sum_part = 0
    for n in orders:
        chi_G_sum_part += 1/n
    
    chi_G = chi_G_sum_part - (num_groups - 1)

    print("\nThe Euler characteristic of G is calculated as: chi(G) = (1/|G1| + 1/|G2| + ...) - (num_groups - 1).")
    print(f"chi(G) = (1/{n1} + 1/{n2} + 1/{n3}) - ({num_groups} - 1) = {chi_G_sum_part} - {num_groups - 1} = {chi_G}")

    # Step 4: Use the formula chi(K) = [G:K] * chi(G).
    # Since K is a free group of rank r, its Euler characteristic is chi(K) = 1 - r.
    chi_K = index * chi_G
    
    print("\nUsing the formula chi(K) = [G:K] * chi(G):")
    # 1 - rank = index * chi_G
    # We print each number in the final equation.
    print(f"1 - rank = {index} * {chi_G}")
    print(f"1 - rank = {chi_K}")

    # Step 5: Solve for the rank.
    rank = 1 - chi_K
    
    print("\nSolving for the rank:")
    print(f"rank = 1 - ({chi_K})")
    print(f"rank = {int(rank)}")
    
    return int(rank)

final_rank = solve_rank_of_kernel()
# The final answer is an integer.
# No need to print again as the function already prints the steps.