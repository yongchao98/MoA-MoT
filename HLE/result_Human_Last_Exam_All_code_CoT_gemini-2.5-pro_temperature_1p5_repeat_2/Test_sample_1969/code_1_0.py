import sys

def solve():
    """
    Calculates the Shapley value for player k in the given coalitional game.
    """
    try:
        n_str = input("Enter the total number of people (n > 1): ")
        n = int(n_str)
        if n <= 1:
            print("Error: n must be greater than 1.")
            return

        k_str = input(f"Enter the index of the person k (1 <= k <= {n}): ")
        k = int(k_str)
        if not (1 <= k <= n):
            print(f"Error: k must be between 1 and {n}.")
            return

    except ValueError:
        print("Invalid input. Please enter integers for n and k.")
        return

    # C1 is the sum of the first n integers
    C1 = n * (n + 1) // 2
    
    # C2 is the sum of the first n squares
    C2 = n * (n + 1) * (2 * n + 1) // 6

    # Sum of indices for all players except k
    C1_minus_k = C1 - k
    
    # Sum of squares of indices for all players except k
    C2_minus_k = C2 - k**2

    # Calculate the four terms of the formula
    term1 = (k**2) * (C1**2)
    term2 = k * (C1_minus_k**3)
    term3 = k * C1_minus_k * C2_minus_k
    term4 = 2 * (k**2) * C2_minus_k

    # Calculate the final Shapley value for player k
    c_k = term1 + term2 + term3 + term4

    # Print the breakdown of the formula
    print("\nThe formula for c_k is: k^2*C1^2 + k*(C1-k)^3 + k*(C1-k)*(C2-k^2) + 2*k^2*(C2-k^2)")
    print(f"For n={n} and k={k}:")
    print(f"C1 (sum 1..n) = {C1}")
    print(f"C2 (sum 1^2..n^2) = {C2}")
    print("\nCalculating the terms:")
    print(f"Term 1: {k}^2 * {C1}^2 = {term1}")
    print(f"Term 2: {k} * ({C1}-{k})^3 = {k} * {C1_minus_k}^3 = {term2}")
    print(f"Term 3: {k} * ({C1}-{k}) * ({C2}-{k}^2) = {k} * {C1_minus_k} * {C2_minus_k} = {term3}")
    print(f"Term 4: 2 * {k}^2 * ({C2}-{k}^2) = 2 * {k**2} * {C2_minus_k} = {term4}")

    print("\nThe final equation is:")
    print(f"c_{k} = {term1} + {term2} + {term3} + {term4}")
    print(f"The fair share for person p_{k} (c_{k}) is: {c_k}")

solve()