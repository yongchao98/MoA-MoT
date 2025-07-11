import networkx as nx
from collections import Counter

def check_if_path_graph(graph):
    """
    Checks if a given graph is a path graph.
    A graph is a path if it is connected and for n > 1 nodes,
    two nodes have degree 1 and the rest (n-2) have degree 2.
    A single-node or two-node connected graph is also a path.
    """
    num_nodes = graph.number_of_nodes()
    if num_nodes == 0:
        return False
    # Path graphs must be connected
    if not nx.is_connected(graph):
        return False
    if num_nodes <= 2:
        return True # A single node or two connected nodes form a path

    degrees = [d for _, d in graph.degree()]
    degree_counts = Counter(degrees)

    # For a path graph with >2 nodes, expect two nodes of degree 1 (endpoints)
    # and n-2 nodes of degree 2 (internal nodes).
    return degree_counts.get(1, 0) == 2 and degree_counts.get(2, 0) == num_nodes - 2

def solve_markov_chain_problem():
    """
    Solves the Markov chain conditioning problem by analyzing the graph structure of the PDF.
    """
    print("Step 1: Define the graphical model for the probability distribution.")
    print("An edge exists between two variables if they are coupled in the PDF.")

    nodes = ['x1', 'x2', 'x3', 'x4', 'x5']
    # Edges derived from the factors of the PDF:
    # {x1,x2,x3}, {x3,x4}, {x2,x3,x4}, {x1,x2,x5}
    edges = {
        ('x1', 'x2'), ('x1', 'x3'), ('x2', 'x3'),
        ('x3', 'x4'),
        ('x2', 'x4'),
        ('x1', 'x5'), ('x2', 'x5')
    }

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(list(edges))

    print(f"\nThe full graph has nodes {G.nodes()} and edges {sorted([sorted(e) for e in G.edges()])}.")
    print("\nStep 2: Check each variable as a conditioning variable.")
    print("Conditioning on a variable means removing it from the graph.")
    print("The remaining graph must form a single connected path (a Markov chain).\n")

    solution = None
    
    for var in sorted(nodes):
        print(f"--- Analyzing conditioning on {var} ---")
        
        # Create the conditional graph by removing the node
        conditional_graph = G.copy()
        conditional_graph.remove_node(var)
        
        # Check if the remaining graph is a connected path
        is_markov_chain = check_if_path_graph(conditional_graph)
        
        degrees = sorted(d for _, d in conditional_graph.degree())
        print(f"Resulting graph has nodes {sorted(conditional_graph.nodes())} and degrees {degrees}.")

        if is_markov_chain:
            print(f"Conclusion: Conditioning on {var} results in a valid Markov chain.")
            solution = var
        else:
            print(f"Conclusion: Conditioning on {var} does not result in a valid Markov chain.")
        print("-" * 45)
    
    print("\nStep 3: Final Answer")
    if solution:
        print(f"The only variable that satisfies the condition is {solution}.")
    else:
        print("No single variable satisfies the condition.")

solve_markov_chain_problem()