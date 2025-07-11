import math

def calculate_shapley_value(n, k):
    """
    Calculates the Shapley value for player p_k in the given game.

    Args:
        n (int): The total number of players.
        k (int): The index of the player for whom to calculate the value.
    """
    if not isinstance(n, int) or not isinstance(k, int) or n <= 1 or k < 1 or k > n:
        print("Invalid input. 'n' must be an integer greater than 1, and 'k' must be an integer between 1 and n.")
        return

    # Formula: c_k = k * ( (n * (n+1) / 2) ^ 3 )
    
    # Step 1: Calculate the sum of the first n integers, S = n*(n+1)/2
    sum_n = n * (n + 1) // 2
    
    # Step 2: Cube the sum, S^3
    sum_n_cubed = sum_n ** 3
    
    # Step 3: Multiply by the player index k
    c_k = k * sum_n_cubed

    # Output the explanation and the step-by-step calculation
    print(f"To find the fair division for player p_{k} out of {n} players:")
    print("The general formula is: c_k = k * ( (n * (n+1) / 2) ^ 3 )")
    print("\n--- Calculation Steps ---")
    print(f"1. Calculate the sum of indices from 1 to n:")
    print(f"   Sum = {n} * ({n} + 1) / 2 = {sum_n}")
    
    print(f"\n2. Cube the sum:")
    print(f"   Sum^3 = {sum_n} ^ 3 = {sum_n_cubed}")
    
    print(f"\n3. Multiply by the player's index k = {k}:")
    print(f"   c_{k} = {k} * {sum_n_cubed} = {c_k}")
    
    print("\n--- Final Answer ---")
    print(f"The exact amount of money that player p_{k} gets is: ${c_k}")

# --- User Input ---
# You can change these values for n and k to see the result for a different case.
total_players_n = 5
player_index_k = 4

# Execute the calculation
calculate_shapley_value(total_players_n, player_index_k)