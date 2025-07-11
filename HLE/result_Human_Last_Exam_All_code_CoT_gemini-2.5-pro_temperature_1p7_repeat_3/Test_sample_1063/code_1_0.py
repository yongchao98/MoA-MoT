import collections

def analyze_actor_graph():
    """
    This function defines and analyzes a graph of actors based on specific
    co-starring criteria to determine its properties.
    """
    # Step 1: Define the nodes (the six actors)
    nodes = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    
    # Step 2: Define the edges based on verified filmography data.
    # An edge exists if two actors were in the same TV series/miniseries season
    # that premiered between 2017-2022.
    # - Aaron Ashmore & Emilia Jones: 'Locke & Key' (premiered 2020)
    # - Krysten Ritter & Charlie Cox: 'The Defenders' (premiered 2017)
    # - Devery Jacobs & Thomas Elms: 'The Order' (premiered 2019)
    edges = [
        ("Aaron Ashmore", "Emilia Jones"),
        ("Krysten Ritter", "Charlie Cox"),
        ("Devery Jacobs", "Thomas Elms")
    ]

    # Step 3: Create an adjacency list to represent the graph
    adj_list = {node: [] for node in nodes}
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)

    print("--- Graph Definition ---")
    print(f"The graph has {len(nodes)} nodes (actors) and {len(edges)} edges (relationships).")
    print("The edges are:")
    for u, v in edges:
        print(f"  - An edge between '{u}' and '{v}'")
    print("-" * 25)

    # Step 4: Analyze the graph's properties
    visited = set()
    
    # --- Check for Connectivity ---
    # We perform a single traversal (like BFS) starting from an arbitrary node.
    # If not all nodes are visited, the graph is disconnected.
    q = collections.deque([nodes[0]])
    visited_for_connectivity = {nodes[0]}
    while q:
        u = q.popleft()
        for v in adj_list[u]:
            if v not in visited_for_connectivity:
                visited_for_connectivity.add(v)
                q.append(v)

    is_connected = len(visited_for_connectivity) == len(nodes)
    
    # --- Check for Cycles ---
    # We use a helper function to perform DFS and detect back edges.
    def has_cycle_util(u, parent, visited_for_cycles):
        visited_for_cycles.add(u)
        for v in adj_list[u]:
            if v not in visited_for_cycles:
                if has_cycle_util(v, u, visited_for_cycles):
                    return True
            elif v != parent:
                return True
        return False

    is_cyclic = False
    visited_for_cycles = set()
    for node in nodes:
        if node not in visited_for_cycles:
            if has_cycle_util(node, None, visited_for_cycles):
                is_cyclic = True
                break

    print("--- Graph Analysis Results ---")
    # Output the mathematical properties used in the final answer
    print(f"Number of Nodes = {len(nodes)}")
    print(f"Number of Edges = {len(edges)}")
    print(f"Is the graph connected? {'Yes' if is_connected else 'No'}")
    print(f"Does the graph contain a cycle? {'Yes' if is_cyclic else 'No'}")
    print("-" * 25)

    # Step 5: Determine the final answer based on the analysis
    if not is_connected and not is_cyclic:
        final_answer = "A. Disconnected and acyclic"
    elif not is_connected and is_cyclic:
        final_answer = "B. Disconnected and cyclic"
    elif is_connected and not is_cyclic:
        final_answer = "C. Connected and acyclic"
    elif is_connected and is_cyclic:
        final_answer = "D. Connected and cyclic, but not a cycle graph"
    else: # This assumes the graph would be a cycle graph
        final_answer = "E. A cycle graph"
        
    print(f"Conclusion: The graph is best described as '{final_answer}'.")

if __name__ == '__main__':
    analyze_actor_graph()