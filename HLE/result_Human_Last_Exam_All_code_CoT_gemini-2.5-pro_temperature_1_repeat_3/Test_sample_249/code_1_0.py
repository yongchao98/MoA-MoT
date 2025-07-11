import math

def solve_diameter():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    # Let n and m be positive integers.
    # We use example values here. You can change them to any positive integers.
    n = 8
    m = 4

    # The problem is valid for m >= 2 and n+1 >= m (with the exception of n=0, m=2)
    # Our formula works even for these edge cases.
    if n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return

    print(f"Given n = {n}, m = {m}")
    print(f"The tree has {n+2} vertices and {m} leaves.")
    print("-" * 20)

    # The total "length" to be distributed among m branches is n+1.
    total_length = n + 1
    print(f"Total length to distribute in branches = n + 1 = {total_length}")

    # We distribute this length as evenly as possible.
    # q is the base length of each branch.
    # r is the number of branches that get an extra +1 length.
    q = total_length // m
    r = total_length % m

    print(f"q = floor((n + 1) / m) = floor({total_length} / {m}) = {q}")
    print(f"r = (n + 1) % m = {total_length} % {m} = {r}")
    print("-" * 20)

    # Determine the diameter based on the value of r.
    if r == 0:
        # All branches have length q. Diameter = q + q.
        diameter = 2 * q
        print("Since r = 0, the two longest branches both have length q.")
        print(f"Minimum Diameter = 2 * q = 2 * {q} = {diameter}")
    elif r == 1:
        # One branch has length q+1, the rest have length q. Diameter = (q+1) + q.
        diameter = 2 * q + 1
        print("Since r = 1, the longest branch has length q+1 and the second-longest has length q.")
        print(f"Minimum Diameter = 2 * q + 1 = 2 * {q} + 1 = {diameter}")
    else: # r >= 2
        # At least two branches have length q+1. Diameter = (q+1) + (q+1).
        diameter = 2 * q + 2
        print("Since r >= 2, the two longest branches both have length q+1.")
        print(f"Minimum Diameter = 2 * q + 2 = 2 * {q} + 2 = {diameter}")

solve_diameter()