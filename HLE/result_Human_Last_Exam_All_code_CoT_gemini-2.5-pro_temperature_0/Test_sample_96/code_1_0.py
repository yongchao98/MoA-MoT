import collections

def solve():
    """
    Calculates the number of torsion elements of order 10 and minimal length
    in the group A/Z for the Artin group A of type E8.
    """
    # The E8 Dynkin diagram can be represented by its adjacency list.
    # Nodes are labeled 1 to 8.
    #   2
    #   |
    # 1-3-4-5-6-7-8
    e8_edges = [
        (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)
    ]
    
    adj = collections.defaultdict(list)
    for u, v in e8_edges:
        adj[u].append(v)
        adj[v].append(u)

    nodes = list(range(1, 9))

    # Step 1: Find all subgraphs of type A4 (paths with 4 nodes).
    # We can find them by exploring all simple paths of length 3.
    a4_subgraphs = []
    for start_node in nodes:
        # A path is a list of nodes.
        q = collections.deque([[start_node]])
        while q:
            path = q.popleft()
            if len(path) == 4:
                # Use a sorted tuple to represent the subgraph uniquely
                subgraph = tuple(sorted(path))
                if subgraph not in a4_subgraphs:
                    a4_subgraphs.append(subgraph)
                continue

            last_node = path[-1]
            for neighbor in adj[last_node]:
                if neighbor not in path:
                    new_path = list(path)
                    new_path.append(neighbor)
                    q.append(new_path)

    # Number of Coxeter elements in A4 is 2^(4-1) = 8.
    num_coxeter_elements_a4 = 8
    
    total_elements = 0
    contributions = []

    # Step 2 & 3: For each A4, find commuting simple reflections and calculate the number of elements.
    for a4 in sorted(list(a4_subgraphs)):
        a4_nodes = set(a4)
        
        # Find neighbors of the A4 subgraph
        neighbors = set()
        for node in a4_nodes:
            for neighbor in adj[node]:
                if neighbor not in a4_nodes:
                    neighbors.add(neighbor)
        
        # Find commuting nodes (not in A4 and not neighbors)
        commuting_nodes = []
        for node in nodes:
            if node not in a4_nodes and node not in neighbors:
                commuting_nodes.append(node)
        
        num_commuting_nodes = len(commuting_nodes)
        contribution = num_commuting_nodes * num_coxeter_elements_a4
        contributions.append(contribution)
        total_elements += contribution
        
        print(f"For A4 subgraph with nodes {sorted(list(a4_nodes))}:")
        print(f"  Commuting simple reflections (nodes): {sorted(commuting_nodes)}")
        print(f"  Number of elements = {num_commuting_nodes} * {num_coxeter_elements_a4} = {contribution}")
        print("-" * 20)

    # Step 4: Print the final sum.
    equation_parts = [str(c) for c in contributions]
    print(f"Total number of elements = {' + '.join(equation_parts)} = {total_elements}")

solve()