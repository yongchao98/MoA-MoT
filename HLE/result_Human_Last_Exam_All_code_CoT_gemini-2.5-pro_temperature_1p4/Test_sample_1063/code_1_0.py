import networkx as nx
import itertools

def analyze_actor_graph():
    """
    This function builds and analyzes a graph of actors based on their shared TV series appearances.
    """
    # Step 1: Define the filmography data for each actor for TV series/miniseries
    # with first episodes released between 2017-2022.
    # Format: { Actor: [Set of TV shows] }
    filmography = {
        "Aaron Ashmore": {"Locke & Key", "Cardinal"},
        "Krysten Ritter": {"The Defenders"},
        "Emilia Jones": {"Locke & Key"},
        "Charlie Cox": {"The Defenders"},
        "Devery Jacobs": {"Cardinal", "The Order"},
        "Thomas Elms": {"The Order"}
    }

    actors = list(filmography.keys())
    
    # Initialize the graph
    G = nx.Graph()
    G.add_nodes_from(actors)
    
    print("Finding edges (connections) between actors...")
    
    # Step 2: Identify edges by checking all pairs of actors for shared shows.
    # Use itertools.combinations to get all unique pairs of actors.
    for actor1, actor2 in itertools.combinations(actors, 2):
        # Find the intersection of their TV show sets
        shared_shows = filmography[actor1].intersection(filmography[actor2])
        
        if shared_shows:
            # Add an edge to the graph if they have a shared show
            G.add_edge(actor1, actor2)
            print(f"- Edge found: {actor1} and {actor2} both appeared in {', '.join(shared_shows)}.")
            
    # Step 3: Analyze the graph's properties.
    print("\n--- Graph Analysis ---")
    
    # Check for connectivity
    if nx.is_connected(G):
        connectivity = "Connected"
    else:
        connectivity = "Disconnected"
        
    # Check for cycles. A forest is a graph with no cycles (acyclic).
    if nx.is_forest(G):
        cyclicity = "Acyclic"
    else:
        cyclicity = "Cyclic"

    print(f"The graph is {connectivity}.")
    print(f"The graph is {cyclicity}.")
    
    # Step 4: Determine the correct description
    print("\n--- Conclusion ---")
    if connectivity == "Disconnected" and cyclicity == "Acyclic":
        print("The graph is Disconnected and Acyclic, which corresponds to Choice A.")
    elif connectivity == "Disconnected" and cyclicity == "Cyclic":
        print("The graph is Disconnected and Cyclic, which corresponds to Choice B.")
    elif connectivity == "Connected" and cyclicity == "Acyclic":
        print("The graph is Connected and Acyclic, which corresponds to Choice C.")
    elif connectivity == "Connected" and cyclicity == "Cyclic":
        # Further check if it's a cycle graph
        is_cycle_graph = True
        if len(G.nodes) > 2:
            for node in G.nodes:
                if G.degree(node) != 2:
                    is_cycle_graph = False
                    break
        else:
            is_cycle_graph = False
            
        if is_cycle_graph:
            print("The graph is a Cycle Graph, which corresponds to Choice E.")
        else:
            print("The graph is Connected and Cyclic, but not a cycle graph, which corresponds to Choice D.")

if __name__ == "__main__":
    # The networkx library is required. If not installed, uncomment the line below.
    # import os; os.system('pip install networkx')
    analyze_actor_graph()
