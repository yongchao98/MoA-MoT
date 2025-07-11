import networkx as nx

def analyze_actor_graph():
    """
    This function builds and analyzes a graph of actor collaborations
    to determine its properties.
    """
    # Step 1: Define the actors (nodes)
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]

    # Step 2: Define the connections (edges) based on filmography research
    # for TV series/miniseries released between 2017 and 2022.
    connections = [
        {"actors": ("Aaron Ashmore", "Emilia Jones"), "show": "Locke & Key (2020-2022)"},
        {"actors": ("Krysten Ritter", "Charlie Cox"), "show": "The Defenders (2017)"},
        {"actors": ("Devery Jacobs", "Thomas Elms"), "show": "The Order (2019-2020)"}
    ]

    # Step 3: Build the graph
    G = nx.Graph()
    G.add_nodes_from(actors)
    for conn in connections:
        G.add_edge(conn["actors"][0], conn["actors"][1], title=conn["show"])

    print("Analyzing the graph of actor co-starring roles (2017-2022)...")
    print("\nFound the following connections (edges):")
    for u, v, data in G.edges(data=True):
        print(f"- {u} and {v} in '{data['title']}'")

    # Step 4: Analyze the graph's properties
    print("\nGraph Analysis:")
    is_connected = nx.is_connected(G)
    # An acyclic graph is also known as a forest.
    is_acyclic = nx.is_forest(G)
    num_components = nx.number_connected_components(G)

    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    print(f"Is the graph connected? {'Yes' if is_connected else 'No'}")
    if not is_connected:
        print(f"The graph consists of {num_components} separate components.")
    print(f"Is the graph acyclic (contains no cycles)? {'Yes' if is_acyclic else 'No'}")

    # Step 5: Determine the final answer
    print("\nConclusion:")
    if not is_connected and is_acyclic:
        print("The graph is Disconnected and acyclic.")
        final_answer = "A"
    elif not is_connected and not is_acyclic:
        print("The graph is Disconnected and cyclic.")
        final_answer = "B"
    elif is_connected and is_acyclic:
        print("The graph is Connected and acyclic.")
        final_answer = "C"
    elif is_connected and not is_acyclic:
        # A cycle graph requires every node to have a degree of 2.
        degrees = [d for n, d in G.degree()]
        if all(d == 2 for d in degrees):
            print("The graph is a cycle graph.")
            final_answer = "E"
        else:
            print("The graph is Connected and cyclic, but not a cycle graph.")
            final_answer = "D"
    else:
        # Fallback case
        print("Could not determine the graph type.")
        final_answer = "Unknown"
    
    # This print is for the final answer submission
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    analyze_actor_graph()