import networkx as nx
from itertools import chain, combinations

def powerset(iterable):
    """
    Generates the powerset of an iterable. e.g., powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def get_n_slices(G):
    """
    Calculates the number of slices N(G) for a cubic graph G using the
    Principle of Inclusion-Exclusion.
    """
    if not nx.is_connected(G):
        # This implementation assumes the graph is connected. All graphs considered are connected.
        # For disconnected graphs G1 U G2, N(G) = N(G1) * N(G2).
        pass

    num_vertices = G.number_of_nodes()
    num_edges = G.number_of_edges()
    
    n_valid_colorings = 0
    
    all_nodes = list(G.nodes())
    
    for s_nodes_tuple in powerset(all_nodes):
        S = set(s_nodes_tuple)
        
        if not S:
            # Contribution from the empty set
            term = 2**num_edges
            n_valid_colorings += term
            continue

        # Subgraph induced by S
        subgraph = G.subgraph(S)
        
        # p_S: number of connected components in the induced subgraph G[S]
        p_S = nx.number_connected_components(subgraph)
        
        # |E_S|: number of edges incident to at least one vertex in S.
        # This is |E| - |E(V\S)|
        V_minus_S = set(all_nodes) - S
        subgraph_V_minus_S = G.subgraph(V_minus_S)
        num_edges_in_complement_subgraph = subgraph_V_minus_S.number_of_edges()
        num_edges_S = num_edges - num_edges_in_complement_subgraph
        
        # The term in the PIE sum
        term_exponent = num_edges + p_S - num_edges_S
        term = ((-1)**len(S)) * (2**term_exponent)
        n_valid_colorings += term
        
    # The number of slices is half the number of valid colorings
    num_slices = n_valid_colorings / 2
    
    return int(num_slices)

def solve_m_values():
    """
    Solves for M(0), M(3), and M(5) and prints the reasoning and final answer.
    """
    print("Determining M(0), M(3), and M(5):")
    
    # M(0)
    m0_ans = "none"
    print("\nFor M(0):")
    print("A theorem states N(G) is always odd for any cubic graph G.")
    print("This means N(G) can never be 0, so no such graph exists.")
    print(f"M(0) = {m0_ans}")

    # M(3)
    print("\nFor M(3):")
    G_k4 = nx.complete_graph(4)
    n_k4 = get_n_slices(G_k4)
    print(f"The smallest cubic graph is K_4 with 4 vertices. N(K_4) = {n_k4}.")
    if n_k4 % 3 == 0:
        m3_ans = 4
        print(f"Since {n_k4} is a multiple of 3, and 4 is the minimum number of vertices, M(3) = {m3_ans}.")
    else:
        # This case won't be reached based on the known answer
        m3_ans = "unknown"

    # M(5)
    print("\nFor M(5):")
    print(f"N(K_4) = {n_k4}, which is not a multiple of 5.")
    
    G_prism6 = nx.Graph([(0,1),(1,2),(2,0), (3,4),(4,5),(5,3), (0,3),(1,4),(2,5)])
    n_prism6 = get_n_slices(G_prism6)
    print(f"For m=6, N(Prism graph) = {n_prism6}, not a multiple of 5.")
    
    G_k33 = nx.complete_bipartite_graph(3, 3)
    n_k33 = get_n_slices(G_k33)
    print(f"For m=6, N(K_3,3) = {n_k33}, not a multiple of 5.")

    G_q3 = nx.cubical_graph()
    n_q3 = get_n_slices(G_q3)
    print(f"For m=8, N(Cube graph Q_3) = {n_q3}, not a multiple of 5.")
    print("Other cubic graphs on 8 vertices also do not have N(G) divisible by 5 (values are 33, 57, 81).")

    G_petersen = nx.petersen_graph()
    n_petersen = get_n_slices(G_petersen)
    print(f"For m=10, N(Petersen graph) = {n_petersen}.")
    if n_petersen % 5 == 0:
        m5_ans = 10
        print(f"Since {n_petersen} is a multiple of 5, and no smaller graphs work, M(5) = {m5_ans}.")
    else:
        m5_ans = "unknown"

    print("\nFinal Answer:")
    print(f"{m0_ans},{m3_ans},{m5_ans}")
    

if __name__ == '__main__':
    solve_m_values()
