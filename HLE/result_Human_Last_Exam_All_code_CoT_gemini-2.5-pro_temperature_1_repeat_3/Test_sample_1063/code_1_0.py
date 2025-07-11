import networkx as nx
import matplotlib.pyplot as plt

def solve_graph_problem():
    """
    This function models the actor graph and determines its properties.
    """
    # The six actors are the nodes of our graph.
    nodes = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]

    # Edges are formed if actors co-starred in a TV series/miniseries
    # with a first episode release year between 2017 and 2022.
    # Based on research:
    # 1. Aaron Ashmore & Emilia Jones in 'Locke & Key' (2020)
    # 2. Krysten Ritter & Charlie Cox in 'The Defenders' (2017)
    # 3. Devery Jacobs & Thomas Elms in 'The Order' (2019)
    # 4. Aaron Ashmore & Thomas Elms in 'SkyMed' (2022)
    edges = [
        ("Aaron Ashmore", "Emilia Jones"),
        ("Krysten Ritter", "Charlie Cox"),
        ("Devery Jacobs", "Thomas Elms"),
        ("Aaron Ashmore", "Thomas Elms"),
    ]

    # Create the graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # Analyze the graph properties
    is_connected = nx.is_connected(G)
    
    # A graph is acyclic if it is a forest. A forest is a collection of trees.
    is_acyclic = nx.is_forest(G)

    # Print the findings
    print("Graph Analysis:")
    print(f"Nodes: {list(G.nodes())}")
    print(f"Edges: {list(G.edges())}")
    print("-" * 20)
    print(f"Is the graph connected? {is_connected}")
    print(f"Is the graph acyclic? {is_acyclic}")
    print("-" * 20)

    # Determine the correct description
    if not is_connected and is_acyclic:
        result = "A. Disconnected and acyclic"
    elif not is_connected and not is_acyclic:
        result = "B. Disconnected and cyclic"
    elif is_connected and is_acyclic:
        result = "C. Connected and acyclic"
    elif is_connected and not is_acyclic:
        # Further check for cycle graph
        num_nodes = G.number_of_nodes()
        num_edges = G.number_of_edges()
        is_cycle_graph = True
        if num_nodes > 0 and num_nodes == num_edges:
            for node in G.nodes():
                if G.degree(node) != 2:
                    is_cycle_graph = False
                    break
        else:
            is_cycle_graph = False
            
        if is_cycle_graph:
            result = "E. A cycle graph"
        else:
            result = "D. Connected and cyclic, but not a cycle graph"
    else:
        result = "Unknown graph type"

    print(f"The graph is best described as: {result}")

solve_graph_problem()
# The final answer is A based on the analysis.
print("\n<<<A>>>")