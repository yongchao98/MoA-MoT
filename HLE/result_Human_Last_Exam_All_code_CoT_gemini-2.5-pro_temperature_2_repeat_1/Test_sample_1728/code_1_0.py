def explain_and_calculate(n):
    """
    Explains the logic and calculates the maximal chromatic number
    for a graph G that is the sum of three cycles of length n.
    """
    print(f"--- Analyzing for n = {n} ---")
    if n < 3:
        print("Error: The length of a cycle (n) must be an integer greater than or equal to 3.")
        return

    # Case 1: n is small (n <= 7)
    if n <= 7:
        result = n
        print("Analysis:")
        print(f"  When n is less than or equal to 7, it is possible to construct the graph G to be a complete graph K_{n}.")
        print(f"  To verify this, we can compare the edges required for K_{n} with the maximum edges three C_{n} can provide.")
        k_n_edges = n * (n - 1) // 2
        max_g_edges = 3 * n
        print(f"  Number of edges required for K_{n}: n * (n - 1) / 2 = {k_n_edges}")
        print(f"  Maximum possible edges in G from three C_{n}: 3 * n = {max_g_edges}")
        print(f"  The condition to form K_{n} is that {k_n_edges} <= {max_g_edges}, which simplifies to n - 1 <= 6, or n <= 7. This holds.")
        print(f"  For all n in [3, 7], constructions for K_n exist.")
        print(f"  The chromatic number of a complete graph K_{n} is n.")

        print(f"\nResult for n = {n}:")
        print(f"  Final Equation: maximal_chromatic_number = n")
        print(f"  Calculated Value: {result}")

    # Case 2: n is large (n > 7)
    else:  # n > 7
        result = 6
        print("Analysis:")
        print(f"  When n > 7, G cannot be a complete graph K_{n}. The maximum degree of any vertex in G is at most 6 (2 from each of the three cycles), but a K_{n} would require a degree of n-1, which is greater than 6.")
        print(f"  By Brooks' Theorem, the chromatic number chi(G) is at most its maximum degree, so chi(G) <= 6.")
        print(f"  We check if a chromatic number of 6 is achievable. This can be done by constructing G to contain a K_6 (a 6-clique) as a subgraph.")
        k_6_edges = 6 * (6 - 1) // 2
        print(f"  A K_6 subgraph has {k_6_edges} edges. It can be shown that three cycles of length n can be strategically chosen to create these 15 edges on a specific subset of 6 vertices.")
        print(f"  The existence of a K_6 subgraph means that chi(G) >= chi(K_6) = 6.")
        print(f"  Since the chromatic number is both at most 6 and at least 6, it must be exactly 6.")
        
        print(f"\nResult for n = {n}:")
        print(f"  Final Equation: maximal_chromatic_number = 6")
        print(f"  Calculated Value: {result}")


if __name__ == '__main__':
    print("This program determines the maximal chromatic number of a graph formed by the sum of three cycles of length n.")
    print("The result depends on the value of n.\n")
    
    # Example for the n <= 7 case
    explain_and_calculate(5)
    print("\n" + "="*40 + "\n")

    # Example for the n <= 7 boundary case
    explain_and_calculate(7)
    print("\n" + "="*40 + "\n")
    
    # Example for the n > 7 case
    explain_and_calculate(10)
    
    # Final conclusion based on the analysis
    print("\n\n<<<The maximal chromatic number is n for 3 <= n <= 7, and 6 for n > 7.>>>")
