def calculate_player_share(n, k):
    """
    Calculates the fair share for player p_k in the described coalitional game.

    The formula for the share c_k is derived using the Shapley value and is:
    c_k = k * ( (n * (n+1)) / 2 )^3

    Args:
        n (int): The total number of players (must be > 1).
        k (int): The index of the player (must be between 1 and n).
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: 'n' must be an integer greater than 1.")
        return
    if not isinstance(k, int) or not (1 <= k <= n):
        print(f"Error: 'k' must be an integer between 1 and {n}.")
        return

    # Calculate the sum of the first n integers
    sum_of_indices = n * (n + 1) // 2

    # The value of the exponent in the formula
    exponent = 3

    # Calculate the final share for player p_k
    share_k = k * (sum_of_indices ** exponent)

    # Print the final equation with the numbers plugged in
    print(f"For n={n} and k={k}, the formula for player p_{k}'s share is:")
    print(f"c_{k} = {k} * (({n} * ({n} + 1)) / 2) ^ {exponent}")

    # Print the final calculated amount
    print(f"\nThe exact amount of money that player p_{k} gets is: ${share_k}")


if __name__ == '__main__':
    # --- User-configurable values ---
    # Set the total number of players 'n'
    total_players_n = 5

    # Set the player index 'k' for whom to calculate the share
    player_index_k = 3
    # --- End of configuration ---

    calculate_player_share(total_players_n, player_index_k)
