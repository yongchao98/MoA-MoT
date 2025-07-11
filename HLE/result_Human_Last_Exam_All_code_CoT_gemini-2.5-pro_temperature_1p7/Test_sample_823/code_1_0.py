import collections

def create_grid_graph(k):
    """Creates a k x k grid graph."""
    adj = collections.defaultdict(list)
    edges = set()
    for r in range(k):
        for c in range(k):
            # Edge to the right
            if c + 1 < k:
                adj[(r, c)].append((r, c + 1))
                adj[(r, c + 1)].append((r, c))
                # Use frozenset for undirected edges to avoid duplicates
                edges.add(frozenset([(r, c), (r, c + 1)]))
            # Edge below
            if r + 1 < k:
                adj[(r, c)].append((r + 1, c))
                adj[(r + 1, c)].append((r, c))
                edges.add(frozenset([(r, c), (r + 1, c)]))
    return adj, list(edges)

def greedy_induced_matching(nodes, edges):
    """
    Finds an induced matching using a greedy algorithm.
    The graph is represented by a set of nodes and a list of edges.
    """
    adj = collections.defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
        
    remaining_nodes = set(nodes)
    remaining_edges = set(frozenset(e) for e in edges)
    
    induced_matching = []

    while remaining_edges:
        # Pick an arbitrary edge
        u, v = next(iter(remaining_edges))
        induced_matching.append((u, v))
        
        # Identify vertices to remove: u, v, and all their neighbors
        neighbors_of_u = set(adj[u])
        neighbors_of_v = set(adj[v])
        
        vertices_to_remove = {u, v} | neighbors_of_u | neighbors_of_v
        
        # Update remaining nodes
        remaining_nodes -= vertices_to_remove
        
        # Update remaining edges: remove any edge incident to a removed vertex
        edges_to_remove = set()
        for e in remaining_edges:
            n1, n2 = e
            if n1 in vertices_to_remove or n2 in vertices_to_remove:
                edges_to_remove.add(e)
        remaining_edges -= edges_to_remove

    return induced_matching

def main():
    k = 10
    n = k * k  # Number of vertices
    d = 4      # Maximum degree of an infinite grid

    adj_list, edge_list = create_grid_graph(k)
    nodes = list(adj_list.keys())

    matching = greedy_induced_matching(nodes, edge_list)
    
    # The lower bound from the proof
    lower_bound = n / (2 * d + 2)
    
    print(f"Graph: {k}x{k} Grid")
    print(f"Number of vertices (n): {n}")
    print(f"Maximum degree (d): {d}")
    print("\n--- Theoretical Lower Bound for Induced Matching Size ---")
    print("Formula: n / (2 * d + 2)")
    print(f"Calculation: {n} / (2 * {d} + 2)")
    print(f"Lower bound: {lower_bound}")

    print("\n--- Greedy Algorithm Result ---")
    print(f"Size of found induced matching: {len(matching)}")
    print(f"Is the result consistent with the theory? {len(matching) >= lower_bound}")
    # print("Found matching edges (sample):", matching[:5])

if __name__ == "__main__":
    main()
