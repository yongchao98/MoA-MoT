import math

def calculate_player_share(n, k):
    """
    Calculates the Shapley value for player p_k in the given coalitional game.

    Args:
        n (int): The total number of players, n > 1.
        k (int): The index of the player, 1 <= k <= n.

    Returns:
        float: The calculated amount c_k for player p_k.
    """
    if not isinstance(n, int) or not isinstance(k, int) or n <= 1 or k < 1 or k > n:
        print("Invalid input. 'n' must be an integer greater than 1, and 'k' must be an integer between 1 and n.")
        return

    # These are the coefficients from the derived formula: c_k = A * [k*B - 6k^2]
    # A = n^2 * (n+1)^2 / 24
    # B = (n+2) * (3n+1)
    
    # Calculate each part of the formula
    n_squared = n**2
    n_plus_1_squared = (n + 1)**2
    
    coeff_A = (n_squared * n_plus_1_squared) / 24
    
    term_n_plus_2 = n + 2
    term_3n_plus_1 = 3 * n + 1
    
    coeff_B = term_n_plus_2 * term_3n_plus_1
    
    k_squared = k**2
    
    # Calculate the final value c_k
    c_k = coeff_A * (k * coeff_B - 6 * k_squared)

    # Output the result showing the formula with computed numbers
    print(f"For n={n} and player p_{k}:")
    print(f"The formula for the player's share c_k is: n^2*(n+1)^2/24 * [k*(n+2)*(3n+1) - 6*k^2]")
    print("\nPlugging in the numbers:")
    # Print the detailed breakdown of the calculation
    print(f"c_{k} = ({n}^2 * ({n}+1)^2 / 24) * [{k}*({n}+2)*(3*{n}+1) - 6*{k}^2]")
    print(f"c_{k} = ({n_squared} * {n_plus_1_squared} / 24) * [{k}*({term_n_plus_2})*({term_3n_plus_1}) - 6*{k_squared}]")
    print(f"c_{k} = {coeff_A:.2f} * [{k}*{coeff_B} - 6*{k_squared}]")
    print(f"c_{k} = {coeff_A:.2f} * [{k*coeff_B} - {6*k_squared}]")
    print(f"c_{k} = {coeff_A:.2f} * [{k*coeff_B - 6*k_squared}]")
    print(f"\nThe final amount for player p_{k} is: ${c_k:,.2f}")

# Example usage with a band of 4 people, finding the share for player p_2
if __name__ == '__main__':
    num_players = 4
    player_index = 2
    calculate_player_share(num_players, player_index)
