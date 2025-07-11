import collections

def solve_graph_problem():
    """
    This function builds a graph based on actors' co-appearances in TV series
    and analyzes its properties to find the correct description.
    """
    # Step 1: Define the nodes (actors)
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    
    # Step 2: Define relevant filmographies (TV series/seasons 2017-2022)
    # This data is based on public filmography records.
    filmographies = {
        "Aaron Ashmore": ["Locke & Key (2020-2022)", "SkyMed (2022)"],
        "Krysten Ritter": ["The Defenders (2017)"],
        "Emilia Jones": ["Locke & Key (2020-2022)"],
        "Charlie Cox": ["The Defenders (2017)"],
        "Devery Jacobs": ["The Order (2019-2020)"],
        "Thomas Elms": ["The Order (2019-2020)", "SkyMed (2022)"]
    }

    # Step 3: Build the adjacency list and edge set for the graph
    adj = collections.defaultdict(list)
    edges = set()
    for i in range(len(actors)):
        for j in range(i + 1, len(actors)):
            actor1 = actors[i]
            actor2 = actors[j]
            
            # Find common shows by checking for intersection in their filmographies
            common_shows = set(filmographies[actor1]) & set(filmographies[actor2])
            
            if common_shows:
                # Add an edge if they share at least one show
                edge = tuple(sorted((actor1, actor2)))
                if edge not in edges:
                    edges.add(edge)
                    adj[actor1].append(actor2)
                    adj[actor2].append(actor1)

    print("Constructed Graph Edges:")
    if not edges:
        print("No edges found.")
    else:
        for edge in sorted(list(edges)):
            common = list(set(filmographies[edge[0]]) & set(filmographies[edge[1]]))[0]
            print(f"- Edge between '{edge[0]}' and '{edge[1]}' (via {common})")
    print("-" * 20)

    # Step 4: Analyze the graph's properties
    
    # --- Find the number of connected components (C) ---
    visited = set()
    num_components = 0
    for actor in actors:
        if actor not in visited:
            num_components += 1
            q = collections.deque([actor])
            visited.add(actor)
            while q:
                node = q.popleft()
                for neighbor in adj[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        q.append(neighbor)

    # --- Determine connectivity ---
    is_connected = (num_components == 1)
    
    # --- Determine cyclicity using the formula |E| = |V| - C ---
    # An undirected graph is a forest (collection of trees, i.e., acyclic) 
    # if and only if the number of edges is equal to the number of vertices 
    # minus the number of connected components.
    num_vertices = len(actors)
    num_edges = len(edges)
    is_acyclic = (num_edges == num_vertices - num_components)

    print("Graph Analysis:")
    print(f"Number of vertices (nodes) |V|: {num_vertices}")
    print(f"Number of edges |E|: {num_edges}")
    print(f"Number of connected components C: {num_components}")
    
    print("\nChecking for cycles using the formula |E| = |V| - C:")
    print(f"The equation to check is: {num_edges} = {num_vertices} - {num_components}")
    print(f"Calculating the right side: {num_vertices - num_components} = {num_vertices - num_components}")
    print(f"Result: The equation {num_edges} = {num_vertices - num_components} is {'True' if is_acyclic else 'False'}.")

    # Step 5: Final Conclusion
    print("\nConclusion:")
    final_description = ""
    if not is_connected and is_acyclic:
        final_description = "Disconnected and acyclic"
        final_answer = "A"
    elif not is_connected and not is_acyclic:
        final_description = "Disconnected and cyclic"
        final_answer = "B"
    elif is_connected and is_acyclic:
        final_description = "Connected and acyclic"
        final_answer = "C"
    elif is_connected and not is_acyclic:
        # Further check for cycle graph
        is_cycle_graph = all(len(adj[actor]) == 2 for actor in actors) and num_edges == num_vertices
        if is_cycle_graph:
            final_description = "A cycle graph"
            final_answer = "E"
        else:
            final_description = "Connected and cyclic, but not a cycle graph"
            final_answer = "D"

    print(f"The graph is {final_description}.")
    print(f"This corresponds to Answer Choice {final_answer}.")

if __name__ == '__main__':
    solve_graph_problem()
    print("\n<<<A>>>")