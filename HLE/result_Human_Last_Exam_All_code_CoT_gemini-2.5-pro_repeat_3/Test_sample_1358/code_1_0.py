import collections

def analyze_graph(nodes, edges):
    """
    Analyzes a graph to check if it's connected and if it's a path.
    A graph is a path if it's connected, has |V|-1 edges, and max degree is 2.
    """
    if not nodes:
        return "No nodes remain.", False

    # Build adjacency list
    adj = collections.defaultdict(list)
    num_edges = 0
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        num_edges += 1

    # Check for non-existent nodes in edge list
    for node in nodes:
        if node not in adj:
            adj[node] = []

    # 1. Check connectivity
    q = collections.deque([nodes[0]])
    visited = {nodes[0]}
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                q.append(v)
    
    is_connected = (len(visited) == len(nodes))
    if not is_connected:
        return "Graph is disconnected.", False

    # 2. Check if it's a path
    degrees = [len(adj[n]) for n in nodes]
    
    is_path_graph = True
    # For a path graph with > 1 nodes, two nodes have degree 1, rest have degree 2.
    if len(nodes) > 1:
        if degrees.count(1) != 2 or degrees.count(2) != len(nodes) - 2:
            is_path_graph = False
    # For a single node graph, it's a path of length 0.
    elif len(nodes) == 1:
         if degrees.count(0) != 1:
            is_path_graph = False
            
    if is_path_graph:
        return "Graph is a connected path (forms a Markov chain).", True
    else:
        return f"Graph is connected but not a path (degrees: {sorted(degrees)}).", False

def main():
    """
    Main function to perform the analysis for each variable.
    """
    variables = ['x1', 'x2', 'x3', 'x4', 'x5']
    
    # Dependencies are derived by treating the conditioned variable as a constant
    # in the log-probability formula and identifying the remaining interaction terms.
    conditional_dependencies = {
        'x1': {('x2', 'x3'), ('x3', 'x4'), ('x2', 'x5')},
        'x2': {('x1', 'x3'), ('x3', 'x4'), ('x1', 'x5')},
        'x3': {('x1', 'x2'), ('x1', 'x5'), ('x2', 'x5')}, # x4 is isolated
        'x4': {('x1', 'x2'), ('x1', 'x3'), ('x2', 'x3'), ('x1', 'x5'), ('x2', 'x5')},
        'x5': {('x1', 'x2'), ('x1', 'x3'), ('x2', 'x3'), ('x2', 'x4'), ('x3', 'x4')}
    }

    print("Analyzing the conditional dependency graph for each variable:\n")
    
    successful_conditions = []
    
    for var in variables:
        remaining_nodes = sorted([v for v in variables if v != var])
        edges = conditional_dependencies[var]
        
        print(f"--- Conditioning on {var} ---")
        print(f"Remaining variables: {remaining_nodes}")
        print(f"Interaction edges: {edges}")
        
        description, is_markov_chain = analyze_graph(remaining_nodes, edges)
        print(f"Result: {description}")
        
        if is_markov_chain:
            successful_conditions.append(var)
        print("-" * (25 + len(var)))
        print()
        
    print("\nConclusion:")
    if len(successful_conditions) == 0:
        print("None of the variables, when conditioned on, result in a Markov chain.")
    elif len(successful_conditions) == 1:
        print(f"Only conditioning on {successful_conditions[0]} results in a Markov chain.")
    else:
        print(f"Conditioning on either {' or '.join(successful_conditions)} results in a Markov chain.")

if __name__ == "__main__":
    main()