import collections

def count_elements():
    """
    This function calculates the number of torsion elements of order 10 and minimal length 5
    in the Coxeter group W(E8).
    """

    # The standard E8 Dynkin diagram is represented by an adjacency list.
    # Nodes are labeled 1 to 8.
    # 1-2-3-4-5-6-7
    #       |
    #       8
    e8_adj = {
        1: [2],
        2: [1, 3],
        3: [2, 4, 8],
        4: [3, 5],
        5: [4, 6],
        6: [5, 7],
        7: [6],
        8: [3]
    }
    
    nodes = list(e8_adj.keys())
    total_count = 0
    
    print("Finding torsion elements of order 10 with minimal length.")
    print("Minimal length is 5. These elements are of the form w = w1 * w2, where:")
    print(" - w2 is a simple reflection (length 1, order 2).")
    print(" - w1 is a Coxeter element of an A4 subsystem (length 4, order 5).")
    print(" - w1 and w2 commute and their supports are disjoint.\n")
    
    # Iterate through each simple reflection s_i as w2
    for i in nodes:
        # Determine the orthogonal complement subgraph for s_i
        perp_nodes = [n for n in nodes if n != i and n not in e8_adj[i]]
        
        perp_adj = collections.defaultdict(list)
        for u in perp_nodes:
            for v in e8_adj[u]:
                if v in perp_nodes:
                    perp_adj[u].append(v)
        
        # Count A4 subgraphs (paths with 4 nodes) in the perp graph.
        # An A4 subgraph is a path P_3.
        num_a4_subgraphs = 0
        
        # We find paths of length 3 (4 nodes).
        # A simple path finding algorithm for small graphs:
        paths = []
        for start_node in perp_nodes:
            # BFS/DFS to find paths of length 3
            q = collections.deque([(start_node, [start_node])])
            while q:
                curr, path = q.popleft()
                if len(path) == 4:
                    # Sort to count each path once
                    paths.append(tuple(sorted(path)))
                    continue
                
                for neighbor in perp_adj[curr]:
                    if neighbor not in path:
                        q.append((neighbor, path + [neighbor]))

        # Since a path can be traversed from two ends, each path is found twice.
        # We count unique sorted tuples of nodes in the paths.
        num_a4_subgraphs = len(set(paths))
        
        # Each A4 subsystem has 8 Coxeter elements of length 4.
        count_for_si = num_a4_subgraphs * 8
        
        print(f"Choice for w2 = s{i}:")
        print(f"  - Orthogonal complement nodes: {sorted(perp_nodes)}")
        print(f"  - Number of A4 subgraphs found: {num_a4_subgraphs}")
        print(f"  - Number of choices for w1: {num_a4_subgraphs} * 8 = {count_for_si}")
        
        total_count += count_for_si
    
    print("\n----------------------------------")
    print(f"Total number of such elements = {total_count}")
    print("----------------------------------")

count_elements()