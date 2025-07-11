import itertools

def count_slices(num_vertices, edges):
    """
    Counts the number of slices for a given cubic graph G.
    A slice is a partition of edges into two classes where each vertex
    is incident to at least one edge of each class.
    
    The method works by checking every possible 2-coloring of the edges.
    A coloring is "valid" if at every vertex, the incident edges are not all
    of the same color. The number of slices is half the number of such
    valid colorings.
    """
    num_edges = len(edges)
    
    adj = {i: [] for i in range(num_vertices)}
    for i, edge in enumerate(edges):
        adj[edge[0]].append(i)
        adj[edge[1]].append(i)

    valid_colorings_count = 0
    # Iterate through all 2^|E| edge colorings (0=class 1, 1=class 2)
    for coloring in itertools.product([0, 1], repeat=num_edges):
        is_valid_coloring = True
        for v in range(num_vertices):
            # For a vertex v, its incident edges can't all have the same color.
            # As G is cubic, each vertex has 3 incident edges.
            edge_indices = adj[v]
            color1 = coloring[edge_indices[0]]
            color2 = coloring[edge_indices[1]]
            color3 = coloring[edge_indices[2]]
            
            if color1 == color2 and color2 == color3:
                is_valid_coloring = False
                break
        
        if is_valid_coloring:
            valid_colorings_count += 1
            
    # The number of slices N(G) is half the number of valid colorings,
    # since swapping the two classes (swapping all 0s and 1s) does
    # not change the partition.
    return valid_colorings_count // 2

def solve():
    """
    Calculates M(0), M(3), and M(5) based on graph properties and slice counting.
    """
    # M(0): We need N(G) = 0. This holds if and only if G has a bridge.
    # The smallest simple cubic graph with a bridge is known to have 10 vertices.
    m0 = 10

    # M(3): We test graphs in increasing order of size.
    # Smallest cubic graph: K4 (m=4)
    k4_vertices = 4
    k4_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = count_slices(k4_vertices, k4_edges)  # This will be 9

    if n_k4 % 3 == 0:
        m3 = 4
    else:
        # Would need to continue searching for m=6, 8, ...
        m3 = "Could not determine with small graphs"


    # M(5): Test graphs in increasing order.
    # For m=4, n_k4 = 9, which is not a multiple of 5.
    # Next size m=6. Let's test the prism graph.
    prism_vertices = 6
    prism_edges = [(0, 1), (1, 2), (2, 0),   # Top triangle
                   (3, 4), (4, 5), (5, 3),   # Bottom triangle
                   (0, 3), (1, 4), (2, 5)]   # Connecting edges
    n_prism = count_slices(prism_vertices, prism_edges) # This will be 15
    
    m5 = "Could not determine with small graphs"
    if n_k4 % 5 != 0:
        if n_prism % 5 == 0:
            m5 = 6

    print(f"{m0},{m3},{m5}")

solve()
<<<10,4,6>>>