import networkx as nx
from collections import Counter

def check_is_path(graph):
    """Checks if a graph with 4 nodes is a path graph."""
    if len(graph.nodes) != 4:
        return False
    # A path must be connected. This also satisfies the problem's condition
    # that no variable is left completely independent.
    if not nx.is_connected(graph):
        return False
    # A path on 4 nodes must be a tree, so no cycles.
    try:
        nx.find_cycle(graph)
        return False # Cycle found
    except nx.NetworkXNoCycle:
        # It's a tree (acyclic and connected). For a tree to be a path,
        # the maximum degree must be 2.
        degrees = [val for (node, val) in graph.degree()]
        if max(degrees) <= 2:
            return True
        else:
            return False

def analyze_markov_chain():
    """
    Analyzes the probability distribution to find the conditioning variable
    that results in a Markov chain.
    """
    variables = ['x1', 'x2', 'x3', 'x4', 'x5']
    
    # Cliques identified from the log-probability function terms
    cliques = [
        {'x1', 'x2', 'x3'},
        {'x3', 'x4'},
        {'x2', 'x3', 'x4'},
        {'x1', 'x2', 'x5'}
    ]

    # Create the graph from the cliques
    G = nx.Graph()
    G.add_nodes_from(variables)
    for clique in cliques:
        nodes_in_clique = list(clique)
        for i in range(len(nodes_in_clique)):
            for j in range(i + 1, len(nodes_in_clique)):
                G.add_edge(nodes_in_clique[i], nodes_in_clique[j])

    print("--- Initial Graph Analysis ---")
    print(f"Nodes: {G.nodes()}")
    print(f"Edges: {G.edges()}")
    print("-" * 30)

    found_solution = False
    # Iterate through each variable to test conditioning on it
    for var_to_condition in variables:
        print(f"\nTesting conditioning on {var_to_condition}:")
        
        # Remaining nodes after conditioning
        remaining_nodes = [v for v in variables if v != var_to_condition]
        
        # Subgraph of the remaining nodes
        subgraph = G.subgraph(remaining_nodes)

        is_path = check_is_path(subgraph)

        degrees = [val for (node, val) in subgraph.degree()]
        degree_counts = Counter(degrees)

        print(f"  Remaining nodes: {list(subgraph.nodes())}")
        print(f"  Edges in subgraph: {list(subgraph.edges())}")
        print(f"  Degrees of nodes: {sorted(degrees)}")

        if is_path:
            print(f"  RESULT: Conditioning on {var_to_condition} turns the distribution into a Markov chain.")
            # Find the endpoints (degree 1) to show the path
            endpoints = [node for node, deg in subgraph.degree() if deg == 1]
            path = nx.shortest_path(subgraph, source=endpoints[0], target=endpoints[1])
            path_str = " - ".join(path)
            print(f"  The resulting Markov chain is: {path_str}")
            found_solution = True
        else:
            if not nx.is_connected(subgraph):
                reason = "the resulting graph is not connected."
            elif len(subgraph.edges()) != len(subgraph.nodes()) - 1:
                reason = f"it has {len(subgraph.edges())} edges instead of the required 3 for a 4-node path."
            else:
                reason = "it contains a cycle or a node with degree > 2."

            print(f"  RESULT: Not a Markov chain because {reason}")
            
if __name__ == '__main__':
    analyze_markov_chain()