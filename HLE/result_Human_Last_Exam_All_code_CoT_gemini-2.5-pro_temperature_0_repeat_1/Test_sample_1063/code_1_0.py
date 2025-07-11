import collections

def analyze_graph_properties():
    """
    This function defines and analyzes the graph of actors based on co-starring roles.
    """
    # The six nodes of the graph
    nodes = {
        "Aaron Ashmore": 0,
        "Krysten Ritter": 1,
        "Emilia Jones": 2,
        "Charlie Cox": 3,
        "Devery Jacobs": 4,
        "Thomas Elms": 5
    }
    node_names = {v: k for k, v in nodes.items()}
    num_nodes = len(nodes)

    # Edges are based on co-starring in a TV series/miniseries season starting 2017-2022.
    # 1. (Aaron Ashmore, Emilia Jones) in "Locke & Key" (S1 2020)
    # 2. (Krysten Ritter, Charlie Cox) in "The Defenders" (2017)
    # 3. (Devery Jacobs, Thomas Elms) in "The Order" (S1 2019)
    # 4. (Aaron Ashmore, Devery Jacobs) in "Cardinal" (S2 2018)
    edges = [
        ("Aaron Ashmore", "Emilia Jones"),
        ("Krysten Ritter", "Charlie Cox"),
        ("Devery Jacobs", "Thomas Elms"),
        ("Aaron Ashmore", "Devery Jacobs")
    ]

    # Build adjacency list
    adj = collections.defaultdict(list)
    for u_name, v_name in edges:
        u, v = nodes[u_name], nodes[v_name]
        adj[u].append(v)
        adj[v].append(u)

    # --- Graph Analysis ---

    # 1. Check for connectivity by counting connected components
    visited = [False] * num_nodes
    num_components = 0
    for i in range(num_nodes):
        if not visited[i]:
            num_components += 1
            q = collections.deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)
    is_connected = (num_components == 1)

    # 2. Check for cycles using DFS
    visited_for_cycle = [False] * num_nodes
    is_cyclic = False
    for i in range(num_nodes):
        if not visited_for_cycle[i]:
            # The `stack` for DFS stores (node, parent)
            stack = [(i, -1)]
            path = set()
            
            while stack:
                u, parent = stack.pop()
                if u in path:
                    is_cyclic = True
                    break
                
                path.add(u)
                visited_for_cycle[u] = True

                for v in adj[u]:
                    if v != parent:
                        stack.append((v, u))
            
            if is_cyclic:
                break

    # --- Print Results ---
    print("Graph Analysis based on TV series co-starring roles (2017-2022):")
    print("\nNodes:")
    for name in nodes:
        print(f"- {name}")

    print("\nEdges:")
    for u_name, v_name in edges:
        print(f"- ({u_name}, {v_name})")

    print("\nProperties:")
    print(f"Number of connected components: {num_components}")
    print(f"Is the graph connected? {'Yes' if is_connected else 'No'}")
    print(f"Does the graph contain a cycle? {'Yes' if is_cyclic else 'No'}")

    print("\nConclusion:")
    if not is_connected and not is_cyclic:
        print("The graph is Disconnected and acyclic.")
    elif not is_connected and is_cyclic:
        print("The graph is Disconnected and cyclic.")
    elif is_connected and not is_cyclic:
        print("The graph is Connected and acyclic (a tree).")
    else: # Connected and cyclic
        print("The graph is Connected and cyclic.")

analyze_graph_properties()