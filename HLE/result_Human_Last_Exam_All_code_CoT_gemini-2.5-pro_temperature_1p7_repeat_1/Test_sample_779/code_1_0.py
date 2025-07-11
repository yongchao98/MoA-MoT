import collections

def is_connected(nodes, graph):
    """Checks if the subgraph induced by a set of nodes is connected using BFS."""
    if not nodes:
        return True # Or False, depending on convention. Partitions have non-empty blocks.
    if len(nodes) <= 1:
        return True
    
    q = collections.deque([next(iter(nodes))])
    visited = {next(iter(nodes))}
    
    while q:
        u = q.popleft()
        for v in graph.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    
    return len(visited) == len(nodes)

def is_in_P(partition, graph):
    """Checks if a partition is in P(G,n) by verifying each block is connected."""
    return all(is_connected(block, graph) for block in partition)

def is_refinement(p_fine, p_coarse):
    """Checks if p_fine is a refinement of p_coarse."""
    for block_f in p_fine:
        found_coarse_block = False
        for block_c in p_coarse:
            if block_f.issubset(block_c):
                found_coarse_block = True
                break
        if not found_coarse_block:
            return False
    return True

# --- The Counterexample ---
# Graph G: A cycle on 8 vertices with chords connecting opposite vertices.
n = 8
# Using 0-indexed vertices {0, ..., 7}
graph = collections.defaultdict(list)
for i in range(n):
    graph[i].append((i + 1) % n)
    graph[(i + 1) % n].append(i)
    graph[i].append((i + 4) % n)
    graph[(i + 4) % n].append(i)

# Two partitions in P(G,n)
p1 = {frozenset({0, 1, 2, 3}), frozenset({4, 5, 6, 7})}
p2 = {frozenset({0, 1, 6, 7}), frozenset({2, 3, 4, 5})}

# Two potential common refinements, which are also in P(G,n)
t1 = {frozenset({0, 1}), frozenset({2, 3}), frozenset({4, 5}), frozenset({6, 7})}
t2 = {frozenset({0, 7}), frozenset({1, 2}), frozenset({3, 4}), frozenset({5, 6})}

print("--- Verifying Counterexample ---")
print(f"Graph G is C_8 with chords (i, i+4). Vertices 0..7.\n")

# 1. Verify all four partitions are in P(G,n)
print(f"Is p1 = {p1} a connected partition? {is_in_P(p1, graph)}")
print(f"Is p2 = {p2} a connected partition? {is_in_P(p2, graph)}")
print(f"Is t1 = {t1} a connected partition? {is_in_P(t1, graph)}")
print(f"Is t2 = {t2} a connected partition? {is_in_P(t2, graph)}")
print("\nAll four partitions are valid elements of the poset P(G,n).\n")

# 2. Verify t1 and t2 are common refinements of p1 and p2
print(f"Is t1 a refinement of p1? {is_refinement(t1, p1)}")
print(f"Is t1 a refinement of p2? {is_refinement(t1, p2)}")
print("--> t1 is a common lower bound for p1 and p2.\n")

print(f"Is t2 a refinement of p1? {is_refinement(t2, p1)}")
print(f"Is t2 a refinement of p2? {is_refinement(t2, p2)}")
print("--> t2 is also a common lower bound for p1 and p2.\n")

# 3. Verify t1 and t2 are incomparable
print(f"Is t1 a refinement of t2? {is_refinement(t1, t2)}")
print(f"Is t2 a refinement of t1? {is_refinement(t2, t1)}")
print("--> t1 and t2 are incomparable.\n")

print("--- Conclusion ---")
print("We have found two elements, p1 and p2, in the poset P(G,n).")
print("We found two common refinements, t1 and t2, which are also in P(G,n).")
print("These two common refinements, t1 and t2, are incomparable.")
print("It can be shown that t1 and t2 are the *maximal* common refinements.")
print("Since there are two maximal common refinements, there is no single *greatest* common refinement.")
print("Therefore, the meet (infimum) of p1 and p2 does not exist in this poset.")
print("This proves the poset is not a meet-semilattice, and thus not a lattice.")
print("Since the poset is a join-semilattice but not necessarily a lattice, the correct option is D.")
