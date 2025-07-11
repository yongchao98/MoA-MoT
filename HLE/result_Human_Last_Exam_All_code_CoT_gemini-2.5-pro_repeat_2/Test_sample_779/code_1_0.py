import collections

def is_connected(graph, nodes):
    """Checks if the subgraph induced by a set of nodes is connected using BFS."""
    if not nodes:
        return True  # An empty set of nodes is trivially connected.
    nodes = set(nodes)
    start_node = next(iter(nodes))
    
    q = collections.deque([start_node])
    visited = {start_node}
    
    while q:
        u = q.popleft()
        for v in graph.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
                
    return visited == nodes

def is_in_P(graph, partition):
    """Checks if a partition is in P(G,n) by checking connectivity of each block."""
    for block in partition:
        if not is_connected(graph, block):
            return False
    return True

def refines(p1, p2):
    """Checks if partition p1 refines partition p2 (p1 <= p2)."""
    p1_blocks = [frozenset(b) for b in p1]
    p2_blocks = [frozenset(b) for b in p2]
    
    for b1 in p1_blocks:
        found_superset = False
        for b2 in p2_blocks:
            if b1.issubset(b2):
                found_superset = True
                break
        if not found_superset:
            return False
    return True

def main():
    # Let n=4 and G be a path graph 1-2-3-4
    n = 4
    # Adjacency list representation of the graph
    G = {
        1: [2],
        2: [1, 3],
        3: [2, 4],
        4: [3]
    }

    # Consider two partitions
    rho1 = [{1, 2}, {3}, {4}]
    rho2 = [{1}, {2, 3}, {4}]

    print(f"Graph G has vertices [1, 2, 3, 4] and edges (1,2), (2,3), (3,4).")
    print(f"Let rho1 = {rho1}")
    print(f"Let rho2 = {rho2}")
    
    # Check if they are in P(G,n)
    is_rho1_in_P = is_in_P(G, rho1)
    is_rho2_in_P = is_in_P(G, rho2)
    
    print(f"\nIs rho1 a connected partition? {is_rho1_in_P}")
    print(f"Is rho2 a connected partition? {is_rho2_in_P}")

    if is_rho1_in_P and is_rho2_in_P:
        # Check for comparability
        rho1_refines_rho2 = refines(rho1, rho2)
        rho2_refines_rho1 = refines(rho2, rho1)

        print(f"\nDoes rho1 refine rho2? {rho1_refines_rho2}")
        print(f"Does rho2 refine rho1? {rho2_refines_rho1}")

        if not rho1_refines_rho2 and not rho2_refines_rho1:
            print("\nConclusion: rho1 and rho2 are incomparable.")
            print("Therefore, the poset is not a total order.")
        else:
            print("\nSomething is wrong in the analysis, rho1 and rho2 should be incomparable.")
    else:
        print("\nSomething is wrong in the analysis, both partitions should be connected.")

if __name__ == "__main__":
    main()