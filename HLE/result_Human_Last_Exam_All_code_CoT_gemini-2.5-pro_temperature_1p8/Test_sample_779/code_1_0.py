import collections

def get_all_partitions(s):
    """Generates all partitions of a set."""
    s = list(s)
    if not s:
        yield []
        return
    first = s[0]
    for smaller_partitions in get_all_partitions(s[1:]):
        for i, subset in enumerate(smaller_partitions):
            yield smaller_partitions[:i] + [subset + [first]] + smaller_partitions[i+1:]
        yield [[first]] + smaller_partitions

def is_connected(nodes, graph):
    """Checks if a set of nodes induces a connected subgraph in the graph."""
    if not nodes:
        return True
    q = collections.deque([list(nodes)[0]])
    visited = {list(nodes)[0]}
    while q:
        u = q.popleft()
        for v in graph.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    return len(visited) == len(nodes)

def to_canonical(p):
    """Converts a partition (list of lists) to a canonical form (tuple of frozensets)."""
    return tuple(sorted([frozenset(b) for b in p], key=lambda fs: min(fs)))

def main():
    n = 4
    vertices = set(range(1, n + 1))
    
    # Define the graph G as a path P_4
    # Edges: 1-2, 2-3, 3-4
    graph_p4 = {
        1: [2],
        2: [1, 3],
        3: [2, 4],
        4: [3]
    }
    
    print("Graph G is a path on 4 vertices: 1-2-3-4.")

    all_partitions_raw = get_all_partitions(vertices)
    all_partitions = {to_canonical(p) for p in all_partitions_raw}

    # P(G,n) is the set of connected partitions
    P_G_n = set()
    for p in all_partitions:
        is_p_connected = all(is_connected(block, graph_p4) for block in p)
        if is_p_connected:
            P_G_n.add(p)
    
    print(f"\nThere are {len(P_G_n)} connected partitions in P(P4, 4).")
    
    # Define the lattice operations
    def leq(p1, p2):
        """Check if p1 <= p2 (p1 is a refinement of p2)."""
        for block1 in p1:
            found_super_block = False
            for block2 in p2:
                if block1.issubset(block2):
                    found_super_block = True
                    break
            if not found_super_block:
                return False
        return True

    def join(p1, p2):
        """Computes the join of two partitions."""
        adj = collections.defaultdict(list)
        nodes = set()
        for p in [p1, p2]:
            for block in p:
                nodes.update(block)
                if len(block) > 1:
                    it = iter(block)
                    u = next(it)
                    for v in it:
                        adj[u].append(v)
                        adj[v].append(u)

        visited = set()
        join_blocks = []
        for i in nodes:
            if i not in visited:
                component = []
                q = collections.deque([i])
                visited.add(i)
                while q:
                    u = q.popleft()
                    component.append(u)
                    for v in adj.get(u, []):
                        if v not in visited:
                            visited.add(v)
                            q.append(v)
                join_blocks.append(frozenset(component))
        return tuple(sorted(join_blocks, key=lambda fs: min(fs)))

    def meet(p1, p2, P, leq_func):
        """Computes the meet of two partitions within the lattice P."""
        lower_bounds = {p for p in P if leq_func(p, p1) and leq_func(p, p2)}
        # The meet is the join of all lower bounds
        if not lower_bounds:
            return None # Should not happen as bot is always there
        
        current_meet = list(lower_bounds)[0]
        for p in list(lower_bounds)[1:]:
            current_meet = join(current_meet, p)
        return current_meet
    
    # A lattice is semimodular if for all x, y: (x covers x^y) => (xVy covers y)
    # We will check this condition.
    
    found_counterexample = False
    
    # rank(p) = n - |p|
    # p1 covers p2 if p2 < p1 and rank(p1) = rank(p2) + 1
    
    sorted_P = sorted(list(P_G_n), key=len, reverse=True)

    for x in sorted_P:
        if found_counterexample: break
        for y in sorted_P:
            if found_counterexample: break
            
            m = meet(x, y, P_G_n, leq)
            j = join(x, y)
            
            # Rank function
            rank_x = n - len(x)
            rank_y = n - len(y)
            rank_m = n - len(m)
            rank_j = n - len(j)
            
            # Check if x covers its meet with y
            if leq(m, x) and rank_x == rank_m + 1:
                # If semimodular, join should cover y
                if not (leq(y,j) and rank_j == rank_y + 1):
                    found_counterexample = True
                    print("\nFound a counterexample to semimodularity:")
                    print(f"Let x = {x}")
                    print(f"Let y = {y}")
                    print(f"Meet (x ^ y) = {m}")
                    print(f"Join (x v y) = {j}")
                    print(f"\nChecking the covering property 'x covers x^y => xVy covers y':")
                    print(f"  rank(x) = {rank_x}, rank(x^y) = {rank_m}. Since rank(x) = rank(x^y) + 1, x covers x^y.")
                    print(f"  rank(y) = {rank_y}, rank(xVy) = {rank_j}. Since rank(xVy) != rank(y) + 1, xVy does not cover y.")
                    print("\nThis violates the condition for a semimodular lattice.")
                    print("Therefore, the lattice is not geometric.")
                    
    if not found_counterexample:
        print("\nCould not find a counterexample. The lattice might be semimodular.")
        
    print("\nBased on this analysis, P(G,n) is a lattice but not necessarily a geometric one.")
    
main()
<<<C>>>