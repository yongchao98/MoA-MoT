import networkx as nx

def analyze_actor_graph():
    """
    This function creates and analyzes a graph of actors based on co-starring roles.
    """
    actors = {
        "AA": "Aaron Ashmore",
        "KR": "Krysten Ritter",
        "EJ": "Emilia Jones",
        "CC": "Charlie Cox",
        "DJ": "Devery Jacobs",
        "TE": "Thomas Elms"
    }

    # Create a graph
    G = nx.Graph()
    G.add_nodes_from(actors.values())

    # Edges based on TV series/seasons from 2017-2022
    # The "equation" here is the set of connections (edges)
    connections = [
        ("Emilia Jones", "Aaron Ashmore", "Locke & Key (2020)"),
        ("Krysten Ritter", "Charlie Cox", "The Defenders (2017)"),
        ("Devery Jacobs", "Thomas Elms", "The Order (2019)"),
        ("Aaron Ashmore", "Devery Jacobs", "Killjoys (2019)")
    ]

    print("Analyzing graph connections...\n")
    for actor1, actor2, series in connections:
        G.add_edge(actor1, actor2)
        print(f"Edge found: {actor1} --- {actor2} (due to {series})")
    
    print("\n--- Graph Analysis ---")

    # Check for connectivity
    is_connected = nx.is_connected(G)
    print(f"Is the graph connected? {is_connected}")

    # Check for cycles
    try:
        # nx.find_cycle returns an edge list for the first cycle found
        cycle = nx.find_cycle(G)
        is_cyclic = True
    except nx.NetworkXNoCycle:
        is_cyclic = False
    
    print(f"Does the graph contain a cycle? {is_cyclic}")
    
    print("\n--- Conclusion ---")
    if not is_connected and not is_cyclic:
        print("The graph is Disconnected and acyclic.")
        final_answer = "A"
    elif not is_connected and is_cyclic:
        print("The graph is Disconnected and cyclic.")
        final_answer = "B"
    elif is_connected and not is_cyclic:
        print("The graph is Connected and acyclic.")
        final_answer = "C"
    elif is_connected and is_cyclic:
        if len(G.nodes) == len(G.edges): # Basic check for a simple cycle graph
             # Check if every node has degree 2 for a true cycle graph
            is_cycle_graph = all(d == 2 for n, d in G.degree())
            if is_cycle_graph:
                print("The graph is a cycle graph.")
                final_answer = "E"
            else:
                print("The graph is Connected and cyclic, but not a cycle graph.")
                final_answer = "D"
        else:
             print("The graph is Connected and cyclic, but not a cycle graph.")
             final_answer = "D"

    print(f"\nThe correct description is: {final_answer}")
    return final_answer

if __name__ == '__main__':
    analyze_actor_graph()
    print("<<<A>>>")
