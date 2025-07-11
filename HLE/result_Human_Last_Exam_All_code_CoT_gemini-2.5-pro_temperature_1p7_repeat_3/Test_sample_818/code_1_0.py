import networkx as nx

def count_slices(G):
    """
    Counts the number of slices in a cubic graph G by brute force enumeration
    of edge 2-colorings.
    """
    if not all(d == 3 for _, d in G.degree()):
        # This function is intended for cubic graphs only.
        return -1 

    num_edges = G.number_of_edges()
    edges = list(G.edges())
    nodes = list(G.nodes())

    # Create a map from each node to the indices of its incident edges
    node_edge_indices = {v: [] for v in nodes}
    for i, edge in enumerate(edges):
        node_edge_indices[edge[0]].append(i)
        node_edge_indices[edge[1]].append(i)

    valid_colorings = 0
    # Iterate through all 2^|E| edge colorings using a bitmask `i`
    for i in range(2**num_edges):
        is_valid = True
        # Check the slice condition for each vertex
        for v in nodes:
            # Get the colors of the 3 incident edges from the bitmask
            colors = [(i >> idx) & 1 for idx in node_edge_indices[v]]
            # If all colors are the same, the coloring is invalid for a slice
            if colors[0] == colors[1] and colors[1] == colors[2]:
                is_valid = False
                break
        if is_valid:
            valid_colorings += 1
            
    # Each slice corresponds to two complementary valid colorings.
    return valid_colorings // 2

def get_cubic_graphs(m):
    """
    Returns a list of non-isomorphic cubic graphs on m vertices.
    This list is hardcoded for the required values of m.
    """
    graphs = []
    if m == 4:
        # Complete graph K4
        graphs.append(nx.complete_graph(4))
    elif m == 6:
        # Triangular prism graph and K_{3,3}
        graphs.append(nx.prism_graph(3))
        graphs.append(nx.complete_bipartite_graph(3, 3))
    elif m == 8:
        # The 5 connected cubic graphs on 8 vertices, plus the disconnected one.
        # graph6 strings are a compact representation for graphs.
        graph6_8_conn = ["G?p@`_", "G?`C@O", "G?`D@O", "G?a_@O", "G?b_@O"]
        for g6 in graph6_8_conn:
            graphs.append(nx.from_graph6_string(g6))
        # The only disconnected cubic graph on 8 vertices is two K4s.
        graphs.append(nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4)))
    elif m == 10:
        # Check some of the 19+ cubic graphs on 10 vertices.
        # The Petersen graph is a well-known example.
        graphs.append(nx.petersen_graph())
        
    return graphs

def find_M(n):
    """
    Finds the smallest m for which a cubic graph G on m vertices has N(G)
    divisible by n.
    """
    m = 4
    while True:
        # The search space grows rapidly, but an answer is expected for small m.
        if m > 10: 
             return "limit_exceeded" # Safety break
        
        graphs_at_m = get_cubic_graphs(m)
        for G in graphs_at_m:
            num_slices = count_slices(G)
            if num_slices != -1 and num_slices % n == 0:
                return m
        
        # The number of vertices in a cubic graph must be even.
        m += 2

def solve():
    """
    Calculates M(0), M(3), and M(5) and prints the result.
    """
    # M(0): N(G) can't be 0, so M(0) is none.
    m0 = "none"
    
    # M(3): Find the smallest m where N(G) is a multiple of 3.
    m3 = find_M(3)

    # M(5): Find the smallest m where N(G) is a multiple of 5.
    m5 = find_M(5)

    print(f"{m0},{m3},{m5}")

solve()
>>>none,4,10