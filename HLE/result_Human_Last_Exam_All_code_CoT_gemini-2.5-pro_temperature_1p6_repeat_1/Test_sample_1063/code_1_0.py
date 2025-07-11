import collections

def analyze_actor_graph():
    """
    This function builds and analyzes a graph based on actor collaborations.
    """
    # Step 1: Define the nodes (actors)
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    
    # Define the collaborations that meet the criteria
    collaborations = {
        "Locke & Key (2020-2022)": ["Aaron Ashmore", "Emilia Jones"],
        "The Defenders (2017)": ["Krysten Ritter", "Charlie Cox"],
        "The Order (2019-2020)": ["Devery Jacobs", "Thomas Elms"]
    }
    
    # Step 2: Build the graph's adjacency list
    adj = collections.defaultdict(list)
    edges = set()

    print("Building the graph based on collaborations...")
    for show, cast in collaborations.items():
        print(f"\nProcessing: {show}")
        # Create an edge for each pair of actors in a show
        for i in range(len(cast)):
            for j in range(i + 1, len(cast)):
                actor1, actor2 = cast[i], cast[j]
                adj[actor1].append(actor2)
                adj[actor2].append(actor1)
                
                # Use a sorted tuple to represent the edge canonically
                edge = tuple(sorted((actor1, actor2)))
                if edge not in edges:
                    edges.add(edge)
                    print(f"- Found edge: {actor1} -- {actor2}")

    print("\n--- Final Graph Representation ---")
    print("Adjacency List:")
    for actor in actors:
        # We check if the actor has any connections, otherwise show an empty list
        connections = adj.get(actor, [])
        print(f"- {actor}: {connections}")

    # Step 3: Analyze the graph's properties
    
    # 3a: Check for Connectivity by finding connected components
    visited = set()
    num_components = 0
    for actor in actors:
        if actor not in visited:
            num_components += 1
            stack = [actor]
            visited.add(actor)
            while stack:
                current_node = stack.pop()
                for neighbor in adj[current_node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
    
    is_connected = (num_components == 1)

    # 3b: Check for Cycles using DFS
    visited_for_cycle = set()
    is_cyclic = False
    for actor in actors:
        if actor not in visited_for_cycle:
            # The 'parent' tracking in this DFS prevents misinterpreting
            # the immediate back-edge in an undirected graph as a cycle.
            stack = [(actor, None)]  # (node, parent)
            while stack:
                current_node, parent = stack.pop()
                if current_node in visited_for_cycle:
                    # If we revisit a node in the same traversal, it's a cycle
                    is_cyclic = True
                    break
                visited_for_cycle.add(current_node)
                
                for neighbor in adj[current_node]:
                    if neighbor != parent:
                        stack.append((neighbor, current_node))
            if is_cyclic:
                break
    
    # Step 4: Output the final analysis
    print("\n--- Graph Analysis Results ---")
    print(f"Number of nodes: {len(actors)}")
    print(f"Number of edges: {len(edges)}")
    print(f"Number of connected components: {num_components}")
    
    prop1 = "Disconnected" if not is_connected else "Connected"
    prop2 = "Acyclic" if not is_cyclic else "Cyclic"
    
    print(f"\nThe graph is {prop1} and {prop2}.")

    # Determine the correct answer choice
    if prop1 == "Disconnected" and prop2 == "Acyclic":
        answer = "A. Disconnected and acyclic"
    elif prop1 == "Disconnected" and prop2 == "Cyclic":
        answer = "B. Disconnected and cyclic"
    elif prop1 == "Connected" and prop2 == "Acyclic":
        answer = "C. Connected and acyclic"
    elif prop1 == "Connected" and prop2 == "Cyclic":
        answer = "D. Connected and cyclic, but not a cycle graph" # or E
    else: # Fallback
        answer = "Unknown"

    print(f"This description matches choice: {answer}")

# Execute the analysis
analyze_actor_graph()