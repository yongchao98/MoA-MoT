import itertools

def get_cubic_graphs(m):
    """
    Provides a list of non-isomorphic cubic graphs for small m.
    Each graph is represented by its number of vertices and a list of edges.
    """
    graphs = []
    if m == 4:
        # K_4 (Complete graph on 4 vertices)
        edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        graphs.append({'m': 4, 'edges': edges})
    elif m == 6:
        # Prism graph (Y_3 or C_3 x K_2)
        edges1 = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
        graphs.append({'m': 6, 'edges': edges1})
        # K_3,3 (Utility graph)
        edges2 = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
        graphs.append({'m': 6, 'edges': edges2})
    # Can be extended for m=8, etc. if needed.
    return graphs

def calculate_N(graph):
    """
    Calculates N(G) for a given graph G.
    N(G) is the number of "slices", which is half the number of (1,2)-spanning subgraphs.
    """
    m = graph['m']
    edges = graph['edges']
    num_edges = len(edges)
    
    count_1_2_subgraphs = 0
    
    # Iterate through all 2^|E| subgraphs using a bitmask
    for i in range(1, 1 << num_edges):
        degrees = [0] * m
        
        for j in range(num_edges):
            if (i >> j) & 1:
                u, v = edges[j]
                degrees[u] += 1
                degrees[v] += 1
        
        # A (1,2)-spanning subgraph must have all vertex degrees in {1, 2}.
        is_1_2_subgraph = all(1 <= deg <= 2 for deg in degrees)
        
        if is_1_2_subgraph:
            count_1_2_subgraphs += 1
            
    # N(G) = k/2 where k is the number of (1,2)-spanning subgraphs.
    # The number of such subgraphs must be even for cubic graphs.
    return count_1_2_subgraphs // 2

def solve():
    """
    Determines M(0), M(3), and M(5) and prints the result.
    """
    # M(0): Based on graph theory, N(G) > 0 for any cubic graph G.
    # Thus, no graph G has N(G) = 0.
    M0 = "none"

    # Find M(3)
    M3 = "none"
    m = 4
    while M3 == "none":
        graphs = get_cubic_graphs(m)
        if not graphs:
            # We have exhausted our library of small graphs
            break
        for g in graphs:
            N_G = calculate_N(g)
            if N_G > 0 and N_G % 3 == 0:
                M3 = m
                break
        if M3 != "none":
            break
        m += 2

    # Find M(5)
    M5 = "none"
    m = 4
    while M5 == "none":
        graphs = get_cubic_graphs(m)
        if not graphs:
            break
        for g in graphs:
            N_G = calculate_N(g)
            if N_G > 0 and N_G % 5 == 0:
                M5 = m
                break
        if M5 != "none":
            break
        m += 2
        
    # Print the final answer in the required format
    print(f"{M0},{M3},{M5}")

solve()