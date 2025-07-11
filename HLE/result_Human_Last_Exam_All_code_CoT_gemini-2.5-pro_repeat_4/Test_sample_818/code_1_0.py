import itertools

def get_edges(adj):
    """Converts an adjacency list to an edge list."""
    edges = []
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            if u < v:
                edges.append((u, v))
    return edges

def count_valid_edge_subsets(num_vertices, edges):
    """Counts the number of edge subsets where every vertex has degree 1 or 2."""
    num_edges = len(edges)
    k = 0
    # Iterate through all 2^|E| possible subgraphs
    for i in range(1 << num_edges):
        subset_edges = []
        for j in range(num_edges):
            if (i >> j) & 1:
                subset_edges.append(edges[j])
        
        degrees = [0] * num_vertices
        for u, v in subset_edges:
            degrees[u] += 1
            degrees[v] += 1
        
        is_valid = True
        if not subset_edges: # Empty set is not a valid partition component
            is_valid = False

        if is_valid:
            for deg in degrees:
                if deg not in [1, 2]:
                    is_valid = False
                    break
        
        if is_valid:
            k += 1
    return k

def solve():
    """
    Determines M(n) for n=0, 3, 5.
    M(n) is the smallest m for which a cubic graph G with m vertices has N(G) a multiple of n.
    """
    
    # M(0): Find smallest m where N(G) is a multiple of 0 (i.e., N(G)=0).
    # A theorem by Laczkovich states that every cubic graph has a slice, so N(G) > 0.
    # Therefore, M(0) does not exist.
    m_0 = "none"

    # M(3): Find smallest m where N(G) is a multiple of 3.
    # Start with the smallest cubic graph, K4 (m=4).
    m_k4 = 4
    adj_k4 = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    edges_k4 = get_edges(adj_k4)
    k_k4 = count_valid_edge_subsets(m_k4, edges_k4)
    n_k4 = k_k4 // 2
    
    m_3 = "none"
    if n_k4 > 0 and n_k4 % 3 == 0:
        m_3 = 4

    # M(5): Find smallest m where N(G) is a multiple of 5.
    m_5 = "none"
    
    # Check m=4
    if n_k4 > 0 and n_k4 % 5 == 0:
        m_5 = 4
    
    # Check m=6 if not found
    if m_5 == "none":
        # Graph 1: Prism graph (C3 x K2)
        m_prism = 6
        # The vertices are 0-1-2 (top) and 3-4-5 (bottom). Edges are (0,3), (1,4), (2,5).
        adj_prism = [[1, 2, 3], [0, 2, 4], [0, 1, 5], [0, 4, 5], [1, 3, 5], [2, 3, 4]]
        edges_prism = get_edges(adj_prism)
        k_prism = count_valid_edge_subsets(m_prism, edges_prism)
        n_prism = k_prism // 2
        if n_prism > 0 and n_prism % 5 == 0:
            m_5 = 6

        # Graph 2: K3,3
        if m_5 == "none":
            m_k33 = 6
            adj_k33 = [[3, 4, 5], [3, 4, 5], [3, 4, 5], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
            edges_k33 = get_edges(adj_k33)
            k_k33 = count_valid_edge_subsets(m_k33, edges_k33)
            n_k33 = k_k33 // 2
            if n_k33 > 0 and n_k33 % 5 == 0:
                m_5 = 6

    # Check m=8 if not found
    if m_5 == "none":
        # There are 5 cubic graphs with 8 vertices. We test them in order.
        graphs_m8 = [
            # 1. Q3, the cube graph
            [[1, 3, 4], [0, 2, 5], [1, 3, 6], [0, 2, 7], [0, 5, 7], [1, 4, 6], [2, 5, 7], [3, 4, 6]],
            # 2. The "prism-like" graph CUB(8,2)
            [[1, 2, 7], [0, 3, 4], [0, 3, 5], [1, 2, 6], [1, 5, 6], [2, 4, 7], [3, 4, 7], [0, 5, 6]],
            # 3. The graph CUB(8,3)
            [[1, 3, 7], [0, 2, 4], [1, 3, 5], [0, 2, 6], [1, 5, 6], [2, 4, 7], [3, 4, 7], [0, 5, 6]],
            # 4. The graph CUB(8,4)
            [[1, 4, 5], [0, 2, 6], [1, 3, 7], [2, 4, 6], [0, 3, 7], [0, 6, 7], [1, 3, 5], [2, 4, 5]],
            # 5. The graph CUB(8,5)
            [[1, 3, 4], [0, 2, 5], [1, 3, 6], [0, 2, 7], [0, 5, 6], [1, 4, 7], [2, 4, 7], [3, 5, 6]]
        ]
        
        for adj in graphs_m8:
            m_g = 8
            edges_g = get_edges(adj)
            k_g = count_valid_edge_subsets(m_g, edges_g)
            n_g = k_g // 2
            if n_g > 0 and n_g % 5 == 0:
                m_5 = 8
                break
    
    print(f"{m_0},{m_3},{m_5}")

solve()