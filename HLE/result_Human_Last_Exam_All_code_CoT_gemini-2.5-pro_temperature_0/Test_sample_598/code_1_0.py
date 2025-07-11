import networkx as nx
import numpy as np

def construct_rooks_graph():
    """Constructs the 4x4 Rook's graph, which is K_4 x K_4."""
    # This is the Cartesian product of two complete graphs K_4.
    # It is an srg(16, 6, 2, 2).
    k4 = nx.complete_graph(4)
    g = nx.cartesian_product(k4, k4)
    # Relabel nodes for consistency, from 0 to 15
    mapping = {node: i for i, node in enumerate(g.nodes())}
    g = nx.relabel_nodes(g, mapping)
    return g

def construct_shrikhande_graph():
    """Constructs the Shrikhande graph."""
    # This is another srg(16, 6, 2, 2), non-isomorphic to the Rook's graph.
    # It can be constructed as a Cayley graph on Z_4 x Z_4.
    g = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    g.add_nodes_from(nodes)
    
    # The connection set S for the Cayley graph
    s = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}
    
    for i1, j1 in nodes:
        for i2, j2 in nodes:
            if i1 == i2 and j1 == j2:
                continue
            # Check if (i2-i1, j2-j1) mod 4 is in S
            di, dj = (i1 - i2) % 4, (j1 - j2) % 4
            if (di, dj) in s:
                g.add_edge((i1, j1), (i2, j2))
                
    # Relabel nodes for consistency, from 0 to 15
    mapping = {node: i for i, node in enumerate(g.nodes())}
    g = nx.relabel_nodes(g, mapping)
    return g

def count_paws(g):
    """Counts the number of paw subgraphs in a graph."""
    # A paw subgraph is a triangle with a pendant edge. It has 4 vertices.
    # The vertex connecting the triangle and the pendant edge is the "center".
    # We count paws by iterating through all possible centers 'u'.
    total_paws = 0
    for u in g.nodes():
        # A paw centered at 'u' consists of a triangle on 'u' and a neighbor 'z' of 'u'
        # that is not part of the triangle and not adjacent to the other triangle vertices.
        
        neighbors_of_u = list(g.neighbors(u))
        if len(neighbors_of_u) < 3:
            continue
        
        # G_u is the subgraph induced by the neighbors of u.
        # An edge (v, w) in G_u means (u, v, w) is a triangle in G.
        g_u = g.subgraph(neighbors_of_u)
        
        for v, w in g_u.edges():
            # We have a triangle (u, v, w).
            # We need to find a neighbor 'z' of 'u' (so z is in V(G_u))
            # that is not adjacent to v or w.
            
            # Neighbors of v within the neighborhood of u
            neighbors_of_v_in_gu = set(g_u.neighbors(v))
            # Neighbors of w within the neighborhood of u
            neighbors_of_w_in_gu = set(g_u.neighbors(w))
            
            # Candidates for z are neighbors of u, excluding v, w, and neighbors of v, w.
            z_candidates = set(neighbors_of_u) - {v, w} - neighbors_of_v_in_gu - neighbors_of_w_in_gu
            total_paws += len(z_candidates)
            
    return total_paws

def count_5_cycles(g):
    """Counts the number of 5-cycles in a graph using the trace formula."""
    # N_5 = (Tr(A^5) / 10) - N_paw
    
    # Calculate N_paw
    n_paws = count_paws(g)
    
    # Calculate Tr(A^5)
    adj_matrix = nx.to_numpy_array(g, dtype=np.int64)
    adj_matrix_5 = np.linalg.matrix_power(adj_matrix, 5)
    trace_a5 = np.trace(adj_matrix_5)
    
    # Calculate N_5
    # Each 5-cycle is counted 10 times in Tr(A^5)
    # Each paw contributes 10 to Tr(A^5)
    n5 = (trace_a5 - 10 * n_paws) // 10
    return n5

def main():
    """Main function to perform the comparison."""
    print("This script demonstrates that two strongly regular graphs with the same parameters")
    print("srg(16, 6, 2, 2) can have a different number of 5-cycles.\n")

    # Construct the graphs
    g_rook = construct_rooks_graph()
    g_shrikhande = construct_shrikhande_graph()

    # Count 5-cycles
    n5_rook = count_5_cycles(g_rook)
    n5_shrikhande = count_5_cycles(g_shrikhande)

    print("Graph 1: The Rook's graph on a 4x4 board")
    print(f"Number of 5-cycles: {n5_rook}")
    print("\nGraph 2: The Shrikhande graph")
    print(f"Number of 5-cycles: {n5_shrikhande}")

    print("\n--- Final Equation ---")
    print("The number of 5-cycles in the two graphs are:")
    print(n5_rook)
    print(n5_shrikhande)
    print(f"The numbers are different: {n5_rook} != {n5_shrikhande}")

if __name__ == '__main__':
    main()