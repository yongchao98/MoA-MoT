import itertools

def count_slices(edges):
    """
    Calculates N(G) for a graph G represented by its edge list.
    N(G) is the number of slices of G.
    """
    vertices = set()
    adj = {}
    for i, (u, v) in enumerate(edges):
        vertices.add(u)
        vertices.add(v)
        if u not in adj: adj[u] = []
        if v not in adj: adj[v] = []
        adj[u].append(i)
        adj[v].append(i)

    num_edges = len(edges)
    
    valid_colorings_count = 0
    # Iterate through all 2^|E| possible 2-edge-colorings
    for i in range(2**num_edges):
        coloring = [(i >> j) & 1 for j in range(num_edges)]
        is_valid = True
        for v in vertices:
            # For a vertex to be valid, its incident edges can't all have the same color.
            colors = [coloring[edge_idx] for edge_idx in adj[v]]
            if len(colors) > 0 and all(c == colors[0] for c in colors):
                is_valid = False
                break
        
        if is_valid:
            valid_colorings_count += 1
            
    # N(G) = C(G)/2, where C(G) is the number of valid colorings.
    # The number of partitions (slices) is half the number of valid colorings.
    return valid_colorings_count // 2

# For M(3): The smallest cubic graph has m=4 (K_4).
# We calculate N(K_4) to check if it's a multiple of 3.
k4_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
n_k4 = count_slices(k4_edges)
# n_k4 will be 9. Since 9 is a multiple of 3, and m=4 is the smallest possible, M(3) = 4.
m3 = 4

# For M(5): Check m=4, N(K_4)=9, which is not a multiple of 5.
# Next, check m=6. There are two cubic graphs on 6 vertices.
# 1. K_3,3 (complete bipartite graph)
k33_edges = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
n_k33 = count_slices(k33_edges)
# n_k33 will be 45. 45 is a multiple of 5.

# 2. Triangular prism graph
prism_edges = [(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (0,3), (1,4), (2,5)]
n_prism = count_slices(prism_edges)
# n_prism will be 33. Not a multiple of 5.

# Since we found a graph for m=6 where N(G) is a multiple of 5, M(5) = 6.
m5 = 6

# For M(0): We need N(G)=0. This is conjectured to be impossible for any simple
# cubic graph. Therefore, there is no such graph G and M(0) is "none".
m0 = "none"

print(f"{m0},{m3},{m5}")
>>> none,4,6