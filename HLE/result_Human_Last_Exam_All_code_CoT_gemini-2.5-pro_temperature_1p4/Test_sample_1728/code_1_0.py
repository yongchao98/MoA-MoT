def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G,
    where G is the sum of three cycles of length n.
    """
    print("The graph G is the sum of three cycles of length n: G = C_n + C_n + C_n.")
    print("The chromatic number of a sum of graphs is the sum of their chromatic numbers.")
    print("Therefore, chi(G) = chi(C_n) + chi(C_n) + chi(C_n) = 3 * chi(C_n).")
    print("-" * 30)

    # Step 1: Define chi(C_n) for even and odd n
    chi_cn_even = 2
    chi_cn_odd = 3
    print("The chromatic number of a cycle C_n depends on the parity of n (for n >= 3):")
    print(f"If n is even, chi(C_n) = {chi_cn_even}.")
    print(f"If n is odd, chi(C_n) = {chi_cn_odd}.")
    print("-" * 30)

    # Step 2: Calculate chi(G) for both cases
    chi_g_even = 3 * chi_cn_even
    chi_g_odd = 3 * chi_cn_odd
    print("Now we calculate chi(G) for both cases:")
    print(f"Case 1: n is even.")
    print(f"chi(G) = 3 * chi(C_n) = 3 * {chi_cn_even} = {chi_g_even}.")
    print(f"\nCase 2: n is odd.")
    print(f"chi(G) = 3 * chi(C_n) = 3 * {chi_cn_odd} = {chi_g_odd}.")
    print("-" * 30)
    
    # Step 3: Determine the maximal chromatic number
    maximal_chi = max(chi_g_even, chi_g_odd)
    print("The maximal chromatic number is the largest possible value for chi(G).")
    print(f"Comparing the two outcomes ({chi_g_even}, {chi_g_odd}), the maximum is {maximal_chi}.")
    print("\nThis maximum is achieved when n is odd.")
    print("The final calculation for the maximal chromatic number is:")
    print(f"{chi_cn_odd} + {chi_cn_odd} + {chi_cn_odd} = {maximal_chi}")

solve_chromatic_number()