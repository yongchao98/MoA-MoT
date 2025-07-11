import networkx as nx
from networkx.algorithms import isomorphism

def check_graph_automorphism(edges, name):
    """Builds a graph and checks the size of its automorphism group."""
    G = nx.Graph()
    G.add_edges_from(edges)
    
    num_edges = G.number_of_edges()
    num_vertices = G.number_of_nodes()
    
    # Graphs must be simple and connected for this problem
    if not nx.is_connected(G):
        print(f"Skipping '{name}': Graph is not connected.")
        return False, num_edges
    
    # Calculate the size of the automorphism group.
    # For NetworkX versions < 3.0, use gm.group_order
    # For NetworkX versions >= 3.0, this is how you get the group
    gm = isomorphism.GraphMatcher(G, G)
    aut_group_size = sum(1 for _ in gm.isomorphisms_iter())

    print(f"Candidate: {name}")
    print(f"Number of vertices (n): {num_vertices}")
    print(f"Number of edges (e): {num_edges}")
    print(f"Size of automorphism group |Aut(γ)|: {aut_group_size}")
    
    if aut_group_size == 3:
        print("\nFound a graph with |Aut(γ)| = 3.")
        return True, num_edges
    else:
        print("This graph does not meet the condition |Aut(γ)| = 3.")
        print("-" * 30)
        return False, num_edges

def find_smallest_graph():
    """
    Tests a sequence of graphs with increasing number of edges to find the smallest
    one with an automorphism group of size 3.
    """
    print("Starting search for the smallest graph γ with |Aut(γ)| = 3.\n")
    print("Based on group theory, the number of edges 'e' is likely a multiple of 3.")
    print("-" * 30)

    # Candidate for e=3: Star graph K_1,3
    # |Aut(K_1,3)| = S_3 = 6
    e_3_edges = [(0, 1), (0, 2), (0, 3)]
    found, e = check_graph_automorphism(e_3_edges, "Star Graph (e=3)")
    if found: return e

    # Candidate for e=6: C3 with a pendant edge on each vertex
    # |Aut(G)| = D_3 = 6
    e_6_edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5)]
    found, e = check_graph_automorphism(e_6_edges, "Decorated C3 (e=6)")
    if found: return e

    # Candidate for e=9: C3 with a path of length 2 on each vertex
    # |Aut(G)| = D_3 = 6
    e_9_edges = [
        (0, 1), (1, 2), (2, 0),  # C3 core
        (0, 3), (3, 4),          # Path on vertex 0
        (1, 5), (5, 6),          # Path on vertex 1
        (2, 7), (7, 8)           # Path on vertex 2
    ]
    found, e = check_graph_automorphism(e_9_edges, "Decorated C3 with P2 (e=9)")
    if found: return e
    
    # Candidate for e=12: The above e=9 graph with additional 'chiral' edges
    # This construction breaks the reflectional symmetry, leaving only C3.
    e_12_edges = e_9_edges + [
        (4, 5), # Actually, this connects (x1, w2) in theoretical model.
                # The nodes are indexed 0..8.
                # v = {0,1,2}, w={3,5,7}, x={4,6,8}
                # (v0-w0-x0) = (0-3-4), (v1-w1-x1)=(1-5-6), (v2-w2-x2)=(2-7-8)
                # To make it chiral, connect x0-w1, x1-w2, x2-w0
                # Edges: (4,5), (6,7), (8,3)
        (4,5), (6,7), (8,3)
    ]
    found, e = check_graph_automorphism(e_12_edges, "Chiral Decorated C3 (e=12)")
    if found: return e
    
    return None

# --- Main Execution ---
smallest_e = find_smallest_graph()
if smallest_e is not None:
    print(f"\nThe smallest number of edges found is {smallest_e}.")
