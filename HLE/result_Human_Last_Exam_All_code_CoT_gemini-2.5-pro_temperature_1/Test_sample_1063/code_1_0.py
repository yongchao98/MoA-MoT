import collections

def solve_graph_problem():
    """
    This function constructs and analyzes a graph of actor collaborations.
    """
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]

    # Collaborations in TV series/miniseries starting 2017-2022
    # Data gathered from public sources (e.g., IMDb, Wikipedia)
    collaborations = {
        "Locke & Key (2020-2022)": ["Aaron Ashmore", "Emilia Jones"],
        "The Defenders (2017)": ["Krysten Ritter", "Charlie Cox"],
        "Jessica Jones, S2 (2018)": ["Krysten Ritter", "Devery Jacobs"],
        "The Order (2019-2020)": ["Devery Jacobs", "Thomas Elms", "Aaron Ashmore"],
    }

    # Step 1: Build the adjacency list for the graph
    adj_list = {actor: set() for actor in actors}
    for show, cast in collaborations.items():
        for i in range(len(cast)):
            for j in range(i + 1, len(cast)):
                actor1 = cast[i]
                actor2 = cast[j]
                adj_list[actor1].add(actor2)
                adj_list[actor2].add(actor1)

    print("Graph Representation (Adjacency List):")
    for actor, neighbors in adj_list.items():
        print(f"- {actor}: {list(neighbors)}")
    print("-" * 20)

    # Step 2: Analyze graph properties

    # Property 1: Connectedness (using BFS)
    def is_connected(graph, all_nodes):
        if not all_nodes:
            return True
        start_node = all_nodes[0]
        queue = collections.deque([start_node])
        visited = {start_node}
        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return len(visited) == len(all_nodes)

    # Property 2: Cyclicity (using DFS)
    def has_cycle_util(graph, node, visited, recursion_stack, parent):
        visited.add(node)
        recursion_stack.add(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                if has_cycle_util(graph, neighbor, visited, recursion_stack, node):
                    return True
            # If neighbor is in recursion stack and is not the parent, it's a back edge -> cycle
            elif neighbor in recursion_stack and neighbor != parent:
                return True
        recursion_stack.remove(node)
        return False

    def is_cyclic(graph, all_nodes):
        visited = set()
        recursion_stack = set()
        for node in all_nodes:
            if node not in visited:
                # Use parent tracking for undirected graphs
                if has_cycle_util(graph, node, visited, recursion_stack, None):
                    return True
        return False

    # Property 3: Is it a "cycle graph"?
    def is_cycle_graph(graph, all_nodes):
        if not is_connected(graph, all_nodes):
            return False
        for node in all_nodes:
            if len(graph[node]) != 2:
                return False
        return True

    connected = is_connected(adj_list, actors)
    cyclic = is_cyclic(adj_list, actors)
    cycle_graph = is_cycle_graph(adj_list, actors)
    
    print("Graph Analysis Results:")
    print(f"Is the graph connected? {connected}")
    print(f"Does the graph have a cycle? {cyclic}")
    print(f"Is the graph a cycle graph? {cycle_graph}")
    print("-" * 20)

    print("Conclusion:")
    if not connected and not cyclic:
        print("The graph is disconnected and acyclic.")
    elif not connected and cyclic:
        print("The graph is disconnected and cyclic.")
    elif connected and not cyclic:
        print("The graph is connected and acyclic.")
    elif connected and cyclic and not cycle_graph:
        print("The graph is connected and cyclic, but not a cycle graph.")
    elif cycle_graph:
        print("The graph is a cycle graph.")

solve_graph_problem()
<<<D>>>