import collections

def analyze_markov_chain_conditioning():
    """
    Analyzes the probability distribution to find which conditioned variable
    results in a Markov chain, following the plan outlined.
    """
    # 1. Define the structure of the graphical model via its maximal cliques
    original_cliques = [
        {'x1', 'x2', 'x3'},
        {'x3', 'x4'},
        {'x1', 'x2', 'x5'}
    ]
    all_vars = {'x1', 'x2', 'x3', 'x4', 'x5'}
    
    conditioning_candidates = sorted(list(all_vars))
    valid_conditions = []

    print("Analyzing the effect of conditioning on each variable...")
    print("-" * 50)

    # 2. Loop through each variable to test conditioning on it
    for cond_var in conditioning_candidates:
        print(f"Testing condition on: {cond_var}")
        
        remaining_vars = all_vars - {cond_var}
        
        # Determine the cliques in the graph of remaining variables
        conditional_cliques = []
        for clique in original_cliques:
            if cond_var in clique:
                new_clique = clique - {cond_var}
                if len(new_clique) >= 2:
                    conditional_cliques.append(new_clique)
            else:
                conditional_cliques.append(clique)

        # Build the conditional graph
        graph = {v: set() for v in remaining_vars}
        for clique in conditional_cliques:
            nodes = list(clique)
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    u, v = nodes[i], nodes[j]
                    graph[u].add(v)
                    graph[v].add(u)
        
        # 3. Check connectivity and if the graph is a path
        is_mc, reason = is_path_graph(graph)
        
        print(f"  Graph of remaining variables {sorted(list(remaining_vars))}:")
        for node, neighbors in sorted(graph.items()):
            print(f"    Node {node} is connected to: {sorted(list(neighbors))}")

        if is_mc:
            print(f"  Result: The conditional graph IS a Markov Chain ({reason}).")
            valid_conditions.append(cond_var)
        else:
            print(f"  Result: The conditional graph IS NOT a Markov Chain ({reason}).")
        print("-" * 50)

    # 5. Conclude based on the results
    print("=== FINAL CONCLUSION ===")
    if 'x1' in valid_conditions and 'x2' in valid_conditions and len(valid_conditions) == 2:
        print("Conditioning on either x1 or x2 results in a connected path graph (a Markov chain).")
        final_answer = 'E'
    elif 'x1' in valid_conditions and len(valid_conditions) == 1:
        final_answer = 'A'
    elif 'x2' in valid_conditions and len(valid_conditions) == 1:
        final_answer = 'B'
    elif 'x3' in valid_conditions and len(valid_conditions) == 1:
        final_answer = 'C'
    elif 'x4' in valid_conditions and len(valid_conditions) == 1:
        final_answer = 'D'
    # Other combinations for completeness
    elif 'x1' in valid_conditions and 'x3' in valid_conditions and len(valid_conditions) == 2:
        final_answer = 'F'
    elif 'x2' in valid_conditions and 'x3' in valid_conditions and len(valid_conditions) == 2:
        final_answer = 'G'
    elif 'x1' in valid_conditions and 'x2' in valid_conditions and 'x3' in valid_conditions and len(valid_conditions) == 3:
        final_answer = 'H'
    else:
        print("No single variable, when conditioned on, satisfies all the criteria.")
        final_answer = 'I'
        
    print(f"\nThe correct option is {final_answer}.")
    return final_answer

def is_path_graph(graph):
    """
    Checks if a graph is a single connected path.
    Returns a tuple (bool, reason_string).
    """
    nodes = list(graph.keys())
    num_nodes = len(nodes)
    
    if num_nodes == 0:
        return True, "Path (empty graph)"

    # Check for connectivity using Breadth-First Search (BFS)
    queue = collections.deque([nodes[0]])
    visited = {nodes[0]}
    while queue:
        u = queue.popleft()
        for v in graph[u]:
            if v not in visited:
                visited.add(v)
                queue.append(v)
    
    if len(visited) != num_nodes:
        return False, "Disconnected graph"

    # A connected graph is a path if it has the correct degree distribution.
    degrees = [len(neighbors) for neighbors in graph.values()]
    if num_nodes == 1: # Single node
        return True, "Path (single node)"
    
    degree_counts = collections.Counter(degrees)
    
    # Path on N>1 nodes must have two nodes of degree 1 and N-2 nodes of degree 2
    if degree_counts.get(1) == 2 and degree_counts.get(2) == num_nodes - 2:
        return True, "Connected with path structure"
    else:
        return False, "Connected, but not a path structure"

if __name__ == '__main__':
    final_answer = analyze_markov_chain_conditioning()
    # The final answer format for the platform
    # print(f"<<<{final_answer}>>>")