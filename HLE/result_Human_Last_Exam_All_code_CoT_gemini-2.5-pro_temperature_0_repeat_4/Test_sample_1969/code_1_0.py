import math

def calculate_fair_share(n, k):
    """
    Calculates the Shapley value for player p_k in the given coalitional game.

    The formula for the Shapley value c_k for player p_k is:
    c_k = k * (sum of all player indices from 1 to n)^(4-1)
    c_k = k * (1 + 2 + ... + n)^3

    Args:
        n (int): The total number of people in the band (n > 1).
        k (int): The index of the person for whom to calculate the share (1 <= k <= n).
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.")
        return
    if not isinstance(k, int) or not (1 <= k <= n):
        print(f"Error: k must be an integer between 1 and {n}.")
        return

    # The formula for c_k is k * (n*(n+1)/2)^3
    print(f"The general formula for the fair share c_k for player p_k is:")
    print("c_k = k * (1 + 2 + ... + n)^3")
    print("Using the formula for the sum of the first n integers, this becomes:")
    print("c_k = k * (n * (n + 1) / 2)^3\n")

    print(f"Calculating the share for player p_{k} with n = {n}:")
    
    # Step 1: Calculate the sum of the first n integers
    sum_of_indices = n * (n + 1) // 2
    print(f"1. Sum of indices from 1 to {n}:")
    print(f"   Sum = {n} * ({n} + 1) / 2 = {sum_of_indices}")

    # Step 2: Cube the sum
    sum_cubed = sum_of_indices ** 3
    print(f"2. Cube the sum:")
    print(f"   ({sum_of_indices})^3 = {sum_cubed}")

    # Step 3: Multiply by k
    c_k = k * sum_cubed
    print(f"3. Multiply by k = {k}:")
    print(f"   c_{k} = {k} * {sum_cubed} = {c_k}")

    print(f"\nThe exact amount of money player p_{k} gets is ${c_k}.")

# Example usage:
# Let's assume there are n=4 people in the band.
# We want to find the share for player p_2 (k=2).
n_example = 4
k_example = 2
calculate_fair_share(n_example, k_example)
