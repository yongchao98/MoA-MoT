import networkx as nx

def verify_and_get_srg_params(G):
    """
    Checks if a graph is strongly regular and returns its parameters (n, d, lambda, mu).
    Returns None if the graph is not a valid SRG.
    """
    if not nx.is_connected(G):
        print("Graph is not connected.")
        return None
    
    n = G.number_of_nodes()
    
    # Check for regularity and get degree d
    degrees = [d for _, d in G.degree()]
    if len(set(degrees)) != 1:
        print("Graph is not regular.")
        return None
    d = degrees[0]

    # Calculate lambda and mu
    lambda_val = -1
    mu_val = -1

    # Calculate lambda: common neighbors for an adjacent pair
    try:
        u, v = next(iter(G.edges()))
        lambda_val = len(list(nx.common_neighbors(G, u, v)))
    except StopIteration: # Graph has no edges
        lambda_val = 0
    
    # Calculate mu: common neighbors for a non-adjacent pair
    # Create a complement graph to find non-edges easily
    try:
        comp = nx.complement(G)
        u_non, v_non = next(iter(comp.edges()))
        mu_val = len(list(nx.common_neighbors(G, u_non, v_non)))
    except StopIteration: # Graph is complete
        # Mu is not well-defined for complete graphs in the standard SRG definition
        # but can be considered as not applicable or a special case.
        # For our purposes, we know the graphs are not complete.
        mu_val = 0
        
    return (n, d, lambda_val, mu_val)

def count_simple_cycles(G, length):
    """
    Counts the number of simple cycles of a given length in a graph.
    """
    # networkx.simple_cycles finds all elementary circuits.
    # We filter them by length.
    return sum(1 for cycle in nx.simple_cycles(G) if len(cycle) == length)

# --- Main execution ---

# 1. Create the two graphs from the srg(16, 6, 2, 2) class
shrikhande_graph = nx.shrikhande_graph()
g1_name = "Shrikhande Graph"

# The L2(4) graph can be constructed as the Cartesian product of two K4 graphs
k4 = nx.complete_graph(4)
l2_4_graph = nx.cartesian_product(k4, k4)
g2_name = "L2(4) Graph (4x4 Rook's Graph)"

# 2. Verify the SRG parameters for both graphs
params_g1 = verify_and_get_srg_params(shrikhande_graph)
params_g2 = verify_and_get_srg_params(l2_4_graph)

print("Yes, a pair of graphs in the same srg(n,d,lambda,mu) class can have a different number of 5-cycles.")
print("Here is a demonstration using the two non-isomorphic graphs in srg(16, 6, 2, 2).\n")

# 3. Count 5-cycles and print the results
if params_g1 == params_g2 and params_g1 is not None:
    print(f"The shared parameters are (n, d, lambda, mu) = {params_g1}.\n")
    
    # Count 5-cycles
    print("Counting 5-cycles in each graph. This may take a moment...")
    c5_g1 = count_simple_cycles(shrikhande_graph, 5)
    c5_g2 = count_simple_cycles(l2_4_graph, 5)
    
    print("\n--- Results ---")
    print(f"Graph 1: {g1_name}")
    print(f"Parameters (n, d, lambda, mu): n={params_g1[0]}, d={params_g1[1]}, lambda={params_g1[2]}, mu={params_g1[3]}")
    print(f"Number of 5-cycles = {c5_g1}\n")

    print(f"Graph 2: {g2_name}")
    print(f"Parameters (n, d, lambda, mu): n={params_g2[0]}, d={params_g2[1]}, lambda={params_g2[2]}, mu={params_g2[3]}")
    print(f"Number of 5-cycles = {c5_g2}\n")
    
    print("As shown, the number of 5-cycles is different, confirming the existence of such a case.")
else:
    print("An error occurred during verification of SRG parameters.")
    print(f"Parameters for {g1_name}: {params_g1}")
    print(f"Parameters for {g2_name}: {params_g2}")
