import itertools

def count_slices(num_vertices, edges):
    """
    Calculates N(G) for a given graph G=(V,E).
    A slice is a partition of edges into two classes where each vertex
    is incident to at least one edge in each class. This is equivalent to
    a 2-edge-coloring where no vertex is monochromatic.
    """
    num_edges = len(edges)
    
    # Create an adjacency list where adj[v] contains indices of edges incident to v.
    adj = [[] for _ in range(num_vertices)]
    for i, edge in enumerate(edges):
        u, v = edge
        adj[u].append(i)
        adj[v].append(i)

    valid_colorings_count = 0
    # Iterate through all 2^|E| possible edge colorings using a bitmask.
    for i in range(1 << num_edges):
        is_valid_slice_coloring = True
        for v_idx in range(num_vertices):
            # Check for a monochromatic vertex.
            if len(adj[v_idx]) > 0:
                first_color = (i >> adj[v_idx][0]) & 1
                is_monochromatic = True
                for edge_idx in adj[v_idx][1:]:
                    if ((i >> edge_idx) & 1) != first_color:
                        is_monochromatic = False
                        break
                if is_monochromatic:
                    is_valid_slice_coloring = False
                    break
        
        if is_valid_slice_coloring:
            valid_colorings_count += 1
            
    # Each slice corresponds to two valid colorings (the coloring and its complement).
    num_slices = valid_colorings_count // 2
    return num_slices

# --- M(0) ---
# N(G)=0 iff G has a bridge. The smallest cubic graph with a bridge has m=10.
# We construct such a graph to verify N(G)=0.
# It's made of two components (H1, H2) of 5 vertices, each with degree sequence (3,3,3,3,2),
# joined by a bridge between their degree-2 vertices.
bridged_graph_10 = {
    "m": 10,
    "edges": [
        # H1 on {0,1,2,3,4}, where vertex 4 has degree 2 in the component
        (0,2), (0,3), (1,2), (1,3), (2,3), (4,0), (4,1),
        # H2 on {5,6,7,8,9}, where vertex 9 has degree 2 in the component
        (5,7), (5,8), (6,7), (6,8), (7,8), (9,5), (9,6),
        # Bridge connecting the two components
        (4,9)
    ]
}
n_bridged = count_slices(bridged_graph_10["m"], bridged_graph_10["edges"])
print(f"For the smallest cubic graph with a bridge (m=10), N(G) = {n_bridged}.")
M0 = 10

# --- M(3) ---
# Smallest cubic graph is K_4 with m=4.
k4_graph = { "m": 4, "edges": [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] }
n_k4 = count_slices(k4_graph["m"], k4_graph["edges"])
print(f"For K_4 (m=4), N(K_4) = {n_k4}.")
M3 = 'none'
if n_k4 % 3 == 0:
    M3 = 4

# --- M(5) ---
M5 = 'none'
# Check m=4
if n_k4 % 5 == 0:
    M5 = 4
    print(f"M(5) found at m=4.")
else:
    # Check m=6
    k33_graph = { "m": 6, "edges": [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)] }
    prism_graph = { "m": 6, "edges": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)] }
    n_k33 = count_slices(k33_graph["m"], k33_graph["edges"])
    n_prism = count_slices(prism_graph["m"], prism_graph["edges"])
    print(f"For K_3,3 (m=6), N(K_3,3) = {n_k33}.")
    print(f"For Prism graph (m=6), N(Y_3) = {n_prism}.")
    if n_k33 % 5 == 0 or n_prism % 5 == 0:
        M5 = 6
        print(f"M(5) found at m=6.")
    else:
        # Check m=8
        cube_graph = { "m": 8, "edges": [(0, 1), (1, 3), (3, 2), (2, 0), (4, 5), (5, 7), (7, 6), (6, 4), (0, 4), (1, 5), (2, 6), (3, 7)] }
        n_cube = count_slices(cube_graph["m"], cube_graph["edges"])
        print(f"For the Cube graph (m=8), N(Q_3) = {n_cube}.")
        if n_cube % 5 == 0:
            M5 = 8

print(f"M(0) = {M0}")
print(f"M(3) = {M3}")
print(f"M(5) = {M5}")
print(f"<<<{M0},{M3},{M5}>>>")