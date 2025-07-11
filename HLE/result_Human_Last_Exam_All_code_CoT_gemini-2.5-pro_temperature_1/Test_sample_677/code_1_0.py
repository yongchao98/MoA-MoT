import collections
from itertools import combinations

def get_super_knight_graph(n, m):
    """Generates the adjacency list for a (3,2) super-knight graph on an n x m board."""
    adj = collections.defaultdict(list)
    moves = [(3, 2), (3, -2), (-3, 2), (-3, -2),
             (2, 3), (2, -3), (-2, 3), (-2, -3)]
    for r in range(n):
        for c in range(m):
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < m:
                    adj[(r, c)].append((nr, nc))
    return adj

def find_k33_subgraph(graph, n, m):
    """
    Finds a K_3,3 subgraph.
    Since the graph is bipartite, we partition nodes by color.
    """
    nodes = list(graph.keys())
    # Partition nodes into two sets based on color (r+c parity)
    set1 = {node for node in nodes if (node[0] + node[1]) % 2 == 0}
    set2 = {node for node in nodes if (node[0] + node[1]) % 2 != 0}

    # Iterate through combinations of 3 nodes in the larger set
    if len(set1) < 3 or len(set2) < 3:
        return None
        
    # We check for 3 nodes in one partition (U) that have 3 common neighbors in the other (V)
    for u_nodes in combinations(set1, 3):
        u1, u2, u3 = u_nodes
        # Find common neighbors
        # An empty graph.get(u,[]) handles disconnected vertices
        neighbors1 = set(graph.get(u1, []))
        neighbors2 = set(graph.get(u2, []))
        neighbors3 = set(graph.get(u3, []))
        
        common_neighbors = list(neighbors1.intersection(neighbors2).intersection(neighbors3))
        
        if len(common_neighbors) >= 3:
            v_nodes = common_neighbors[:3]
            # We found a K_3,3. U={u1,u2,u3}, V={v1,v2,v3}
            return (list(u_nodes), v_nodes)

    for u_nodes in combinations(set2, 3):
        u1, u2, u3 = u_nodes
        neighbors1 = set(graph.get(u1,[]))
        neighbors2 = set(graph.get(u2,[]))
        neighbors3 = set(graph.get(u3,[]))
        
        common_neighbors = list(neighbors1.intersection(neighbors2).intersection(neighbors3))
        
        if len(common_neighbors) >= 3:
            v_nodes = common_neighbors[:3]
            return (list(u_nodes), v_nodes)

    return None

def solve():
    """
    Finds the supremum of nm for planar (3,2)-super-knight graphs with n,m >= 4.
    """
    search_limit = 11
    minimal_non_planar = []

    print("Searching for minimal non-planar boards...")
    for n in range(4, search_limit):
        for m in range(n, search_limit):
            is_supergraph = False
            for min_n, min_m in minimal_non_planar:
                if n >= min_n and m >= min_m:
                    is_supergraph = True
                    break
            if is_supergraph:
                continue

            graph = get_super_knight_graph(n, m)
            k33 = find_k33_subgraph(graph, n, m)
            
            if k33:
                print(f"Found a K_3,3 on a {n}x{m} board. This is a minimal non-planar board.")
                # print(f"  U = {k33[0]}, V = {k33[1]}")
                minimal_non_planar.append((n, m))

    if not minimal_non_planar:
        print(f"No non-planar boards found up to {search_limit-1}x{search_limit-1}.")
        print("The supremum may be infinite.")
        return

    print(f"\nMinimal non-planar boards found: {minimal_non_planar}")

    max_planar_size = 0
    
    # We check boards up to a size slightly larger than the minimal non-planar ones found
    check_limit = search_limit + 5 

    for n in range(4, check_limit):
        for m in range(4, check_limit):
            is_non_planar = False
            for min_n, min_m in minimal_non_planar:
                # Check if (n,m) contains a minimal non-planar board
                if (n >= min_n and m >= min_m) or (n >= min_m and m >= min_n):
                    is_non_planar = True
                    break
            
            if not is_non_planar:
                if n * m > max_planar_size:
                    max_planar_size = n * m

    print(f"\nLargest product nm for a planar board is {max_planar_size}.")
    print("This is the supremum of the set of sizes of planar rectangles.")
    print(f"\nFinal Answer: {max_planar_size}")


solve()