import collections

def solve_n3_parity():
    """
    Calculates the number of centrally symmetric partitions for the 3x3 grid.
    This tells us the parity of P_3.
    """
    n = 3
    nodes = set((i, j) for i in range(1, n + 1) for j in range(1, n + 1))
    
    # Adjacency list for the 3x3 grid graph
    adj = collections.defaultdict(list)
    for r in range(1, n + 1):
        for c in range(1, n + 1):
            if r + 1 <= n: adj[(r, c)].append((r + 1, c))
            if r - 1 >= 1: adj[(r, c)].append((r - 1, c))
            if c + 1 <= n: adj[(r, c)].append((r, c + 1))
            if c - 1 >= 1: adj[(r, c)].append((r, c - 1))

    def is_connected(node_set):
        """Checks if a set of nodes induces a connected subgraph using BFS."""
        if not node_set:
            return True
        q = collections.deque([next(iter(node_set))])
        visited = {next(iter(node_set))}
        while q:
            u = q.popleft()
            for v_neighbor in adj[u]:
                if v_neighbor in node_set and v_neighbor not in visited:
                    visited.add(v_neighbor)
                    q.append(v_neighbor)
        return len(visited) == len(node_set)

    def rotate_180(p):
        """Rotates a vertex 180 degrees around the center of the 3x3 grid."""
        return (n + 1 - p[0], n + 1 - p[1])

    center_node = ((n + 1) // 2, (n + 1) // 2)
    piece_size = (n * n) // 3

    # Find all centrally symmetric connected pieces of size 3.
    # A centrally symmetric piece of odd size must contain the center node.
    symmetric_pieces = []
    other_nodes = sorted(list(nodes - {center_node}))
    
    # A symmetric piece is {center, v, g(v)}
    for i in range(len(other_nodes)):
        v = other_nodes[i]
        gv = rotate_180(v)
        # Avoid double counting pairs
        if v > gv:
            continue
        
        piece = frozenset({center_node, v, gv})
        if len(piece) == piece_size and is_connected(piece):
            symmetric_pieces.append(piece)

    print(f"Found {len(symmetric_pieces)} centrally symmetric pieces for n=3.")
    
    num_symmetric_partitions = 0
    # For each symmetric piece C, check if the remaining graph can be partitioned
    # into a pair of pieces {A, B} such that B is the 180-degree rotation of A.
    for C in symmetric_pieces:
        remaining_nodes = nodes - C
        
        # Partition remaining_nodes into two sets A and B, where B = g(A)
        # We just need to find one connected component A of size 3.
        # Then B is determined, and we check if it's the rest of the nodes.
        
        # Find all possible pieces 'A' in the remaining nodes
        import itertools
        possible_A_sets = itertools.combinations(remaining_nodes, piece_size)
        
        partitions_for_C = 0
        
        # Use a set to store found partitions {A,B} to avoid double counting
        found_partitions = set()

        for a_tuple in possible_A_sets:
            A = frozenset(a_tuple)
            if not is_connected(A):
                continue

            B_expected = frozenset(rotate_180(p) for p in A)
            
            # Check if A and B form a partition of the remaining nodes
            if A.isdisjoint(B_expected) and A.union(B_expected) == remaining_nodes:
                # Store the pair {A, B} in a canonical way (sorted)
                canonical_pair = tuple(sorted([A, B_expected], key=lambda f: sorted(list(f))[0]))
                found_partitions.add(canonical_pair)

        # Each valid partition {A, B} of G-C gives one symmetric partition of G.
        num_symmetric_partitions += len(found_partitions)

    print(f"The number of centrally symmetric partitions for n=3 is {num_symmetric_partitions}.")
    if num_symmetric_partitions % 2 == 0:
        print("This is an even number, so P_3 is even.")
    else:
        print("This is an odd number, so P_3 is odd.")

solve_n3_parity()
