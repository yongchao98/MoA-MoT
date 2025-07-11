import collections

def find_non_cut_points(graph: dict[str, list[str]]) -> list[str]:
    """
    Finds all points in a graph that are not cut points.
    A point is a 'cut point' if its removal disconnects the graph.
    For dendrites, non-cut-points correspond to endpoints.

    Args:
        graph: The graph represented as an adjacency list.

    Returns:
        A list of the non-cut-points (endpoints) of the graph.
    """

    def is_connected_after_removal(start_node: str, excluded_node: str) -> bool:
        """
        Checks if the graph remains connected after a node is removed.
        Uses Breadth-First Search (BFS).
        """
        nodes = set(graph.keys())
        nodes.remove(excluded_node)
        
        if not nodes:
            return True # An empty or single-node graph is connected.

        queue = collections.deque([start_node])
        visited = {start_node}

        while queue:
            current_node = queue.popleft()
            for neighbor in graph.get(current_node, []):
                if neighbor != excluded_node and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        return len(visited) == len(nodes)

    non_cut_points = []
    nodes = list(graph.keys())
    
    if len(nodes) <= 2:
        # In a graph with 0, 1, or 2 nodes, no node can be a cut point.
        return nodes

    for p in nodes:
        # To check if 'p' is a cut point, we remove it and see if the
        # remaining graph is connected. We pick an arbitrary remaining
        # node to start the connectivity search.
        remaining_nodes = [node for node in nodes if node != p]
        start_node = remaining_nodes[0]
        
        # If the graph is still connected after removing p, then p is not a cut point.
        if is_connected_after_removal(start_node, p):
            non_cut_points.append(p)
            
    return non_cut_points

def solve():
    """
    Solves the problem by demonstrating the minimal case.
    """
    print("Based on topological theorems, the number of compactifications is equal to the")
    print("number of endpoints of the remainder space X, which must be a dendrite.")
    print("Any non-degenerate dendrite has at least 2 endpoints.")
    print("\nLet's test the simplest dendrite, an arc, represented as a graph:")
    
    # An arc [0,1] can be represented as a simple path graph, e.g., A-B-C
    arc_graph = {
        'A': ['B'],
        'B': ['A', 'C'],
        'C': ['B']
    }

    endpoints = find_non_cut_points(arc_graph)
    num_endpoints = len(endpoints)
    
    print("Graph model of an arc:", arc_graph)
    print("Endpoints (non-cut-points) found:", endpoints)
    
    # The final "equation" is that the number of compactifications for this X
    # is the number of its endpoints.
    print(f"\nFinal calculation for the minimal case (X = arc):")
    print(f"Number of topologically distinct compactifications = Number of endpoints = {num_endpoints}")

    print("\nSince the minimum possible number of endpoints is 2, and we have found a")
    print("space that achieves this, the smallest number of compactifications is 2.")


if __name__ == '__main__':
    solve()
    print("\n<<<2>>>")