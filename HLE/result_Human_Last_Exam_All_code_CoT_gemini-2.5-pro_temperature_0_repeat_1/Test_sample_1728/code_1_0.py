def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the sum of three cycles of length n.
    """
    print("Step 1: Define the graph G and the 'sum' operation.")
    print("The graph G is the sum of three cycles of length n. In graph theory, the 'sum' of graphs is typically the 'join' operation.")
    print("So, G = C_n + C_n + C_n, where C_n is a cycle of length n.")
    print("-" * 40)

    print("Step 2: Use the formula for the chromatic number of a graph join.")
    print("The chromatic number of a join of graphs is the sum of their individual chromatic numbers.")
    print("Formula: χ(G_1 + G_2 + ... + G_k) = χ(G_1) + χ(G_2) + ... + χ(G_k)")
    print("-" * 40)

    print("Step 3: Apply the formula to graph G.")
    print("For our graph G, this means: χ(G) = χ(C_n) + χ(C_n) + χ(C_n) = 3 * χ(C_n).")
    print("-" * 40)

    print("Step 4: Analyze the chromatic number of a cycle, χ(C_n).")
    print("The value of χ(C_n) for a cycle of length n (where n >= 3) depends on whether n is even or odd:")
    chi_cn_even = 2
    chi_cn_odd = 3
    print(f" - If n is even, C_n is bipartite, so χ(C_n) = {chi_cn_even}.")
    print(f" - If n is odd, C_n is not bipartite and requires 3 colors, so χ(C_n) = {chi_cn_odd}.")
    print("-" * 40)

    print("Step 5: Find the maximal value for χ(G).")
    print("To find the maximal chromatic number of G, we must choose the value of n that maximizes χ(C_n).")
    max_chi_cn = max(chi_cn_even, chi_cn_odd)
    print(f"The maximum value for χ(C_n) is {max_chi_cn}, which occurs when n is an odd number.")
    print("-" * 40)

    print("Step 6: Calculate the final answer.")
    print("We substitute the maximal value of χ(C_n) into our equation for χ(G):")
    factor = 3
    final_result = factor * max_chi_cn
    print(f"Maximal χ(G) = {factor} * {max_chi_cn} = {final_result}")

solve_chromatic_number()