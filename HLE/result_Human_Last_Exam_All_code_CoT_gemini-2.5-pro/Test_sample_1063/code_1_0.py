import collections

def analyze_graph():
    """
    This function models and analyzes the graph of actors based on their collaborations.
    """
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]

    # Collaborations found in TV series/miniseries released between 2017-2022
    collaborations = [
        ("Locke & Key (2020-2022)", ["Aaron Ashmore", "Emilia Jones"]),
        ("The Defenders (2017)", ["Krysten Ritter", "Charlie Cox"]),
        ("The Order (2019-2020)", ["Devery Jacobs", "Thomas Elms"]),
    ]

    # Build the adjacency list for the graph
    adj = collections.defaultdict(list)
    edges = set()
    for series, cast in collaborations:
        # For each pair of actors in the cast, add an edge
        for i in range(len(cast)):
            for j in range(i + 1, len(cast)):
                actor1 = cast[i]
                actor2 = cast[j]
                # Add edge in both directions for an undirected graph
                adj[actor1].append(actor2)
                adj[actor2].append(actor1)
                # Store a canonical representation of the edge to avoid duplicates
                edges.add(tuple(sorted((actor1, actor2))))

    print("--- Graph Analysis ---")
    print("Nodes (Actors):")
    for actor in actors:
        print(f"- {actor}")
    
    print("\nEdges (Collaborations):")
    for actor1, actor2 in edges:
        print(f"- ({actor1}, {actor2})")

    # --- Graph Property Analysis ---

    # 1. Connectivity Check
    visited = set()
    q = collections.deque([actors[0]])
    visited.add(actors[0])
    while q:
        node = q.popleft()
        for neighbor in adj[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)

    is_connected = (len(visited) == len(actors))
    connectivity_status = "Connected" if is_connected else "Disconnected"

    # 2. Cyclicity Check (for undirected graphs)
    visited_for_cycle = set()
    is_cyclic = False

    def has_cycle_util(u, parent):
        nonlocal is_cyclic
        visited_for_cycle.add(u)
        for v in adj[u]:
            if v not in visited_for_cycle:
                if has_cycle_util(v, u):
                    return True
            # If an adjacent vertex is visited and is not parent of current vertex,
            # then there is a cycle
            elif v != parent:
                return True
        return False

    for node in actors:
        if node not in visited_for_cycle:
            if has_cycle_util(node, None):
                is_cyclic = True
                break
    
    cyclicity_status = "Cyclic" if is_cyclic else "Acyclic"

    print(f"\nConnectivity: {connectivity_status}")
    print(f"Cyclicity: {cyclicity_status}")
    
    print("\n--- Conclusion ---")
    if connectivity_status == "Disconnected" and cyclicity_status == "Acyclic":
        print("The graph is Disconnected and Acyclic.")
        final_answer = "A"
    elif connectivity_status == "Disconnected" and cyclicity_status == "Cyclic":
        print("The graph is Disconnected and Cyclic.")
        final_answer = "B"
    elif connectivity_status == "Connected" and cyclicity_status == "Acyclic":
        print("The graph is Connected and Acyclic.")
        final_answer = "C"
    elif connectivity_status == "Connected" and cyclicity_status == "Cyclic":
        # Further check could be added here for cycle graph vs. general cyclic
        print("The graph is Connected and Cyclic.")
        final_answer = "D or E" # For this problem, we don't need to distinguish further
    else:
        final_answer = "Unknown"

    # Return the letter corresponding to the final answer
    return final_answer


if __name__ == "__main__":
    answer = analyze_graph()
    print(f"\nFinal Answer choice is: {answer}")
    print(f"<<<{answer}>>>")
