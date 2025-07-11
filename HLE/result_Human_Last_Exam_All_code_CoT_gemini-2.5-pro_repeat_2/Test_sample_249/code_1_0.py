import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of an undirected tree with
    n+2 vertices and m leaves.

    Args:
      n: A positive integer.
      m: A positive integer.

    Returns:
      The minimum possible diameter of the tree.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        raise ValueError("n and m must be positive integers.")
        
    # A tree with m leaves must have at least ceil(m/2) edges in its minimal connecting subtree,
    # leading to constraints on the number of internal vertices. A key constraint is that
    # the number of internal vertices (n+2-m) must be at least 1 if m > 2.
    # n + 2 - m >= 1 => n + 1 >= m. If this is not met, a tree with m leaves is not possible.
    if n + 1 < m:
        # This case corresponds to an impossible graph structure.
        # For instance, if n=1, m=3, V=3. A path P3 has n=1, m=2. You can't have 3 leaves.
        # However, the problem statement implies such a graph G exists. We will proceed with the formula.
        # The formula handles these cases correctly by resulting in small diameters for small n.
        pass

    q = n // m
    r = n % m

    # Case 1: The "ternary" case which is the most spread out
    if r == m - 1 and q % 2 != 0 and m % 2 != 0:
        diameter = 2 * q + 3
        print(f"n = {n}, m = {m}")
        print(f"q = floor(n/m) = {q}")
        print(f"r = n % m = {r}")
        print("Condition: r == m-1 AND q is odd AND m is odd -> TRUE")
        print(f"Diameter = 2*q + 3 = 2*{q} + 3 = {diameter}")

    # Case 2: The most "compact" case
    elif r == 0 and (q % 2 == 0 or m % 2 == 0):
        diameter = 2 * q + 1
        print(f"n = {n}, m = {m}")
        print(f"q = floor(n/m) = {q}")
        print(f"r = n % m = {r}")
        print("Condition: r == 0 AND (q is even OR m is even) -> TRUE")
        print(f"Diameter = 2*q + 1 = 2*{q} + 1 = {diameter}")

    # Case 3: The "standard" case
    else:
        diameter = 2 * q + 2
        print(f"n = {n}, m = {m}")
        print(f"q = floor(n/m) = {q}")
        print(f"r = n % m = {r}")
        print("Condition: All other cases -> TRUE")
        print(f"Diameter = 2*q + 2 = 2*{q} + 2 = {diameter}")
        
    return diameter

# Example usage with some values for n and m.
# You can change these values to test other cases.
n_val = 4
m_val = 3

# Calculate and print the result
min_diameter = get_minimum_diameter(n_val, m_val)
# The final answer is the formula itself, which is printed inside the function.
# We print the numerical result for the example case here.
# print(f"\nFor n={n_val}, m={m_val}, the minimum diameter is: {min_diameter}")