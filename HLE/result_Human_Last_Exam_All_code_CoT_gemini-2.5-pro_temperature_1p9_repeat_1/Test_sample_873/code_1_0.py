import math

def calculate_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H, where G is K_m with a C5
    removed, and H is K_n with a C4 removed.
    """
    # Let's use example values for m and n as they are not provided in the prompt.
    # The logic holds for any natural numbers m >= 5 and n >= 4.
    m = 8
    n = 10

    if m < 5:
        print("Error: m must be a natural number greater than or equal to 5.")
        return
    if n < 4:
        print("Error: n must be a natural number greater than or equal to 4.")
        return
        
    # Step 1: Calculate Theta(H)
    # H = K_n - C_4 is a perfect graph, so Theta(H) = alpha(H) = 2.
    theta_h = 2

    # Step 2: Calculate Theta(G)
    # G = K_m - C_5 = C_5 + K_{m-5}.
    # For G' = C_5, q = m-5, Theta(G' + K_q) = max(Theta(G'), q)
    theta_c5 = math.sqrt(5)
    theta_g = max(theta_c5, m - 5)

    # Step 3: Calculate the final capacity
    # Theta(G⊠H) = Theta(G) * Theta(H)
    final_capacity = theta_g * theta_h

    print(f"For m = {m} and n = {n}:")
    print(f"The Shannon capacity of H is Θ(H) = {theta_h}.")
    print(f"The Shannon capacity of G is Θ(G) = max(√5, m-5) = max({theta_c5:.4f}, {m-5}) = {theta_g:.4f}.")
    print(f"The total Shannon capacity is Θ(G⊠H) = Θ(G) * Θ(H).")
    print(f"Result = {theta_g:.4f} * {theta_h} = {final_capacity:.4f}")

calculate_shannon_capacity()