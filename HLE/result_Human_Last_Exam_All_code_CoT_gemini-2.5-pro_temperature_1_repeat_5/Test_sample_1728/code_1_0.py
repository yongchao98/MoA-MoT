def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G which is the sum of three cycles of length n.
    """

    print("Let G be a graph that is the sum (or join) of three cycles of length n, denoted as G = C_n + C_n + C_n.")
    print("The chromatic number of a graph join H_1 + H_2 is the sum of their individual chromatic numbers: chi(H_1 + H_2) = chi(H_1) + chi(H_2).")
    print("This means that the set of colors used for each subgraph must be disjoint.")
    print("Therefore, the chromatic number of G is chi(G) = chi(C_n) + chi(C_n) + chi(C_n) = 3 * chi(C_n).")
    print("-" * 30)

    print("The chromatic number of a cycle graph, chi(C_n), depends on the parity of n (assuming n >= 3).")

    # Case 1: n is even
    chi_cn_even = 2
    print(f"\nCase 1: n is even (n = 4, 6, ...)")
    print(f"A cycle of even length is bipartite, so its chromatic number is {chi_cn_even}.")
    chi_g_even = 3 * chi_cn_even
    print(f"In this case, chi(G) = 3 * chi(C_n) = 3 * {chi_cn_even} = {chi_g_even}.")

    # Case 2: n is odd
    chi_cn_odd = 3
    print(f"\nCase 2: n is odd (n = 3, 5, ...)")
    print(f"A cycle of odd length is not bipartite and requires {chi_cn_odd} colors.")
    chi_g_odd = 3 * chi_cn_odd
    print(f"In this case, chi(G) = 3 * chi(C_n) = 3 * {chi_cn_odd} = {chi_g_odd}.")
    print("-" * 30)

    # Determine the maximal chromatic number
    maximal_chi = max(chi_g_even, chi_g_odd)
    print("\nTo find the maximal chromatic number, we compare the results from the two cases.")
    print(f"The possible chromatic numbers for G are {chi_g_even} (for even n) and {chi_g_odd} (for odd n).")
    print(f"The maximum of these two values is {maximal_chi}.")

    print("\nThe final equation for the maximal chromatic number, which occurs when n is odd, is:")
    print(f"Maximal chi(G) = {3} * {chi_cn_odd} = {maximal_chi}")

solve_chromatic_number()
<<<9>>>