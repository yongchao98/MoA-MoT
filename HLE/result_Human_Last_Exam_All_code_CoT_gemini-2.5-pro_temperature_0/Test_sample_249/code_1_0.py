import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    
    Args:
        n: A positive integer.
        m: A positive integer.
    """
    # A tree with V>=3 vertices must have at least 2 leaves.
    # n is positive, so n+2 >= 3.
    if m < 2:
        print("Invalid input: A tree with 3 or more vertices must have at least 2 leaves.")
        return

    # The number of leaves cannot exceed the number of non-internal nodes in a simple tree structure.
    # The number of internal nodes is at least 1 for a tree with diameter > 1.
    if m > n + 1:
        print(f"Invalid input: For n={n}, the number of leaves m cannot exceed n+1={n+1}.")
        return

    # Case: The tree is a star graph (one central node).
    # This occurs when there is only 1 internal node.
    # Number of internal nodes k = n + 2 - m = 1  => m = n + 1.
    if m == n + 1:
        print(f"For n={n}, m={m}, the tree is a star graph.")
        print("The minimum diameter is 2.")
        return

    # Case: The tree's internal nodes form a path of length 1 (two nodes).
    # Number of internal nodes k = n + 2 - m = 2 => m = n.
    if m == n:
        print(f"For n={n}, m={m}, the internal nodes form a path of length 1.")
        print("The minimum diameter is 3.")
        return

    # General cases for m <= n-1
    # Condition to check if the internal tree can be made into a star.
    if 2 * m >= n + 1:
        print(f"For n={n}, m={m}, the condition 2*m >= n+1 holds.")
        print("The internal nodes can form a star-like structure with diameter 2.")
        print("The minimum diameter is 2 + 2 = 4.")
    else:
        # The number of leaves m is too small to form a star-like internal tree.
        # The diameter of the internal tree is forced to be larger.
        numerator = 2 * (n - m + 1)
        denominator = m
        diam_I = math.ceil(numerator / denominator)
        result = 2 + diam_I
        print(f"For n={n}, m={m}, the condition 2*m < n+1 holds.")
        print(f"The minimum diameter is 2 + ceil(2 * (n - m + 1) / m) = 2 + ceil(2 * ({n} - {m} + 1) / {m}) = 2 + ceil({numerator}/{denominator}) = {result}.")

# --- Example Usage ---
# You can change these values
n_val = 8
m_val = 3

get_minimum_diameter(n_val, m_val)