import sys

def get_fair_division(n, k):
    """
    Calculates the Shapley value c_k for player p_k in the given coalitional game.

    The formula for c_k is derived to be:
    c_k = k * S1^3 - k^2 * S1^2 + k * S1 * S2
    where:
    S1 is the sum of the first n integers.
    S2 is the sum of the first n squared integers.
    """
    if not isinstance(n, int) or not isinstance(k, int) or n <= 1 or k < 1 or k > n:
        print("Error: 'n' must be an integer greater than 1, and 'k' must be an integer between 1 and n.", file=sys.stderr)
        return

    # S1 = 1 + 2 + ... + n
    s1 = n * (n + 1) // 2
    
    # S2 = 1^2 + 2^2 + ... + n^2
    s2 = n * (n + 1) * (2 * n + 1) // 6

    # The formula for c_k
    term1 = k * (s1 ** 3)
    term2 = (k ** 2) * (s1 ** 2)
    term3 = k * s1 * s2
    
    result = term1 - term2 + term3

    print(f"For n = {n} and k = {k}:")
    print("The formula for c_k is: k * (S1)^3 - k^2 * (S1)^2 + k * S1 * S2")
    print(f"Where S1 (sum of 1 to n) = {s1}")
    print(f"and S2 (sum of squares from 1 to n) = {s2}")
    print("\nCalculating the terms:")
    print(f"Term 1: {k} * ({s1})^3 = {term1}")
    print(f"Term 2: {k**2} * ({s1})^2 = {term2}")
    print(f"Term 3: {k} * {s1} * {s2} = {term3}")

    print("\nFinal calculation:")
    print(f"c_{k} = {term1} - {term2} + {term3} = {result}")


# Example usage:
# You can change these values to see the results for different n and k
n_val = 4
k_val = 2
get_fair_division(n_val, k_val)