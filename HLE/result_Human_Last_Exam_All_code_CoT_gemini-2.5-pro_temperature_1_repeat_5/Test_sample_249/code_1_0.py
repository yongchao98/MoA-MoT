import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.

    Args:
      n: A positive integer.
      m: A positive integer.

    Returns:
      The minimum possible diameter.
    """
    # A tree must have at least 2 leaves if it has more than 2 vertices.
    # The number of vertices is n+2. If n>=1, V>=3.
    # Also, m must be less than the number of vertices n+2.
    if m < 2 and n > 0:
        # A tree with >2 vertices has >= 2 leaves. n>0 -> V>=3.
        # This case is not possible for a tree, but we handle it.
        return -1 # Or raise an error
    if m >= n + 2:
        # Number of leaves cannot be equal to or greater than the number of vertices
        return -1 # Or raise an error

    # I is the number of internal nodes
    I = n - m + 2

    if I <= 0:
      # This happens if m >= n+2, which is not a valid tree.
      # All vertices are leaves, which is only possible for V=1 or V=2.
      # If V=2 (n=0, m=2), diameter is 1. Our I=0 case. Let's adjust logic.
      # Path P_k has n=k-2 vertices and m=2 leaves. Dia is k-1.
      if n == 0 and m == 2: # Graph is an edge
        print("The graph is a single edge (2 vertices, 2 leaves).")
        print("Diameter = 1")
        return 1
      else:
        # Not a valid simple tree
        return -1

    if I == 1:
        # This occurs when m = n + 1.
        # The internal nodes form a single-vertex tree T_I. D_T_I = 0.
        # This corresponds to a star graph.
        print(f"For n={n}, m={m}: The number of internal nodes is 1 (since m = n+1).")
        print("The structure is a star graph, so the minimum diameter is 2.")
        return 2
    elif I == 2:
        # This occurs when m = n.
        # T_I is an edge. D_T_I = 1.
        # This requires at least 2 leaves in G, so m >= 2, thus n >= 2.
        if n < 2:
            print(f"For n={n}, m={m}: This configuration is not possible for a simple tree.")
            return -1
        print(f"For n={n}, m={m}: The number of internal nodes is 2 (since m = n).")
        print("The internal nodes form a path of length 1, so the minimum diameter is 1 + 2 = 3.")
        return 3
    else: # I >= 3
        # Check if a star-like structure for T_I is possible.
        # This requires m >= number of leaves in T_I, which is I - 1 for a star.
        if 2 * m >= n + 1:
            print(f"For n={n}, m={m}: There are {I} internal nodes and enough leaves (2m >= n+1) to make them a star-like structure.")
            print("The diameter of the internal tree can be 2, so the minimum diameter is 2 + 2 = 4.")
            return 4
        else:
            # Not enough leaves for a star T_I, so diameter is larger.
            # D_G = 2 * ceil((n-m)/m) + 2
            k_num = n - m
            k_den = m
            k = math.ceil(k_num / k_den)
            diameter = 2 * k + 2
            print(f"For n={n}, m={m}: There are {I} internal nodes, but not enough leaves for a compact structure.")
            print(f"The minimum diameter is given by the formula 2 * ceil((n-m)/m) + 2.")
            print(f"Diameter = 2 * ceil(({n}-{m})/{m}) + 2 = 2 * {k} + 2 = {diameter}")
            return diameter

# --- Example Usage ---
# You can change these values to test
n_val = 5
m_val = 2

print(f"Calculating minimum diameter for a tree with n={n_val}, m={m_val}:")
print("Total vertices = n+2 =", n_val + 2)
print("Number of leaves = m =", m_val)
print("-" * 20)
min_diameter = get_minimum_diameter(n_val, m_val)
print("-" * 20)
if min_diameter != -1:
    print(f"Final Answer: The minimum possible diameter is {min_diameter}")
