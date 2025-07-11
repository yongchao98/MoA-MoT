import math

def solve_graph_problem(d):
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.

    Args:
        d (int): An even integer representing the degree of vertex v1.
    """
    if not isinstance(d, int) or d < 0 or d % 2 != 0:
        print("Error: d must be a non-negative even integer.")
        return

    # Step 1: Calculate the total number of edges from the removed vertices {v1, v2, v3}
    # to the rest of the graph G'.
    # deg(v1) = d, deg(v2) = d + 1, deg(v3) = d + 1
    # Total edges = d + (d + 1) + (d + 1)
    total_edges_from_V0 = 3 * d + 2

    # Step 2: Determine the maximum possible number of leaf blocks (L) in G'.
    # Each leaf block requires at least 2 edges from {v1, v2, v3}.
    # L_max = total_edges_from_V0 / 2
    max_leaf_blocks = total_edges_from_V0 / 2

    # Step 3: Calculate the minimal number of edges to add to make G' 2-edge-connected.
    # This is ceil(L_max / 2).
    min_new_edges = math.ceil(max_leaf_blocks / 2)

    # Output the steps of the calculation
    print(f"Given d = {d}:")
    print(f"The degrees of the three vertices are {d}, {d+1}, and {d+1}.")
    print(f"The total number of edges connected to these vertices is {d} + {d+1} + {d+1} = {total_edges_from_V0}.")
    print(f"The maximum number of leaf blocks (L) in G' is {total_edges_from_V0} / 2 = {int(max_leaf_blocks)}.")
    print(f"The minimal number of new edges to make G' 2-edge-connected is ceil(L / 2).")
    print(f"Result = ceil({int(max_leaf_blocks)} / 2) = {min_new_edges}.")


# Example usage with d = 10 (an even integer)
# You can change this value to any other even integer.
d_example = 10
solve_graph_problem(d_example)