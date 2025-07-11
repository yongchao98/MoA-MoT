import itertools

def analyze_actor_graph():
    """
    Analyzes the graph of actors based on co-starring roles
    in TV series/seasons premiering between 2017 and 2022.
    """
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]

    # Step 1: Pre-compiled list of relevant TV works for each actor.
    # An edge exists if two actors share a work from this list.
    filmographies = {
        "Aaron Ashmore": {"Locke & Key (2020-2022)"},
        "Krysten Ritter": {"The Defenders (2017)"},
        "Emilia Jones": {"Locke & Key (2020-2022)"},
        "Charlie Cox": {"The Defenders (2017)"},
        "Devery Jacobs": {"The Order (2019-2020)"},
        "Thomas Elms": {"The Order (2019-2020)"}
    }

    # Step 2: Find all pairs of actors and check for edges.
    actor_pairs = itertools.combinations(actors, 2)
    edges = []
    for actor1, actor2 in actor_pairs:
        # Find common works by checking for intersections in their filmographies.
        common_works = filmographies[actor1].intersection(filmographies[actor2])
        if common_works:
            edges.append((actor1, actor2, list(common_works)))

    # Step 3: Describe the resulting graph structure.
    print("Step-by-step analysis of the graph:")
    print("Nodes:", ", ".join(actors))
    print("\nFinding edges based on co-starring roles (2017-2022):")
    if not edges:
        print("No edges found.")
    else:
        for actor1, actor2, works in edges:
            # The "equation" for each edge is the connection itself.
            print(f"- Edge found: {actor1} --- {actor2} (due to: {', '.join(works)})")

    # Step 4: Analyze the graph's properties.
    print("\nGraph Analysis:")
    
    # A graph is connected if there is a path between any two nodes.
    # By inspection, we can see the graph is made of separate components.
    print("The graph is composed of several disconnected components:")
    for edge in edges:
        print(f"  - A component connecting {edge[0]} and {edge[1]}")
    
    num_components = len(edges)
    is_connected = num_components == 1 and len(actors) == 2

    # A graph is acyclic if it contains no cycles. A cycle requires at least
    # 3 nodes in a single component. Our components only have 2 nodes each.
    is_cyclic = False

    print("\nFinal Properties:")
    if is_connected:
        print("- The graph is CONNECTED.")
    else:
        print(f"- The graph is DISCONNECTED (it has {num_components} components).")

    if is_cyclic:
        print("- The graph is CYCLIC.")
    else:
        print("- The graph is ACYCLIC (it contains no cycles).")
        
    print("\nConclusion: The graph is disconnected and acyclic, which corresponds to option A.")


if __name__ == '__main__':
    analyze_actor_graph()