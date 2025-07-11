def solve_maximal_chromatic_number():
    """
    Calculates and explains the maximal chromatic number of a graph
    that is the sum of three cycles of length n.
    """
    print("This script calculates the maximal chromatic number of a graph G, where G is the sum of three cycles of length n.")
    print("-" * 20)

    # Step 1: Explain the formula for the chromatic number of the graph sum
    print("Step 1: The Chromatic Number of a Graph Sum")
    print("The 'sum' of graphs (G1 + G2) refers to the graph join operation.")
    print("A key property is that the chromatic number of a join of graphs is the sum of their individual chromatic numbers:")
    print("χ(G1 + G2) = χ(G1) + χ(G2)")
    print("For our graph G = C_n + C_n + C_n, this means:")
    print("χ(G) = χ(C_n) + χ(C_n) + χ(C_n) = 3 * χ(C_n)")
    print("-" * 20)

    # Step 2: Explain the chromatic number of a cycle
    print("Step 2: The Chromatic Number of a Cycle (C_n)")
    print("The chromatic number of a cycle C_n depends on its length n (for n ≥ 3):")
    chi_cn_even = 2
    chi_cn_odd = 3
    print(f"- If n is even, χ(C_n) = {chi_cn_even} (since it's a bipartite graph).")
    print(f"- If n is odd, χ(C_n) = {chi_cn_odd} (since it requires a third color).")
    print("-" * 20)

    # Step 3: Maximize the chromatic number
    print("Step 3: Finding the Maximal Chromatic Number")
    print("To find the maximal chromatic number of G, we must use the maximal possible value for χ(C_n).")
    max_chi_cn = max(chi_cn_even, chi_cn_odd)
    print(f"Comparing the two cases, the maximum value for χ(C_n) is {max_chi_cn}.")
    print("-" * 20)

    # Step 4: Final calculation
    print("Step 4: Final Calculation")
    print("We substitute this maximal value back into our formula for χ(G):")
    print("Maximal χ(G) = 3 * max(χ(C_n))")
    final_result = 3 * max_chi_cn
    print("\nThe final equation with the numbers is:")
    print(f"{max_chi_cn} + {max_chi_cn} + {max_chi_cn} = {final_result}")


# Execute the function to print the solution steps and the result.
solve_maximal_chromatic_number()

print("\n<<<9>>>")