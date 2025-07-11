import collections

def get_connected_components(nodes, edges):
    """Finds the connected components of a graph subset."""
    if not nodes:
        return []
    
    adj = collections.defaultdict(list)
    node_set = set(nodes)
    for u, v in edges:
        if u in node_set and v in node_set:
            adj[u].append(v)
            adj[v].append(u)
            
    components = []
    visited = set()
    for node in nodes:
        if node not in visited:
            comp = set()
            q = collections.deque([node])
            visited.add(node)
            while q:
                u = q.popleft()
                comp.add(u)
                for v in adj.get(u, []):
                    if v not in visited:
                        visited.add(v)
                        q.append(v)
            components.append(frozenset(comp))
    return components

def compute_join(p1, p2, n):
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
            parent[root_j] = root_i

    for p in [p1, p2]:
        for block in p:
            if not block: continue
            first_el = next(iter(block))
            for el in block:
                union(first_el, el)
    
    components = collections.defaultdict(set)
    for i in range(1, n + 1):
        components[find(i)].add(i)
    
    return [frozenset(v) for v in components.values()]

def compute_meet(p1, p2, all_edges):
    """Computes the meet of two partitions in P(G,n)."""
    # 1. Compute meet in the full partition lattice Pi_n
    pi_inf_blocks = []
    for b1 in p1:
        for b2 in p2:
            intersection = b1.intersection(b2)
            if intersection:
                pi_inf_blocks.append(intersection)
    
    # 2. Refine by splitting non-G-connected blocks
    meet_partition = []
    for block in pi_inf_blocks:
        components = get_connected_components(list(block), all_edges)
        meet_partition.extend(components)
    return meet_partition

def format_partition(p):
    """Helper function to print partitions nicely."""
    # Sort blocks by their minimum element for consistent ordering
    sorted_blocks = sorted([sorted(list(b)) for b in p])
    return '{' + ', '.join(['{' + ', '.join(map(str, b)) + '}' for b in sorted_blocks]) + '}'

# Define the graph G and n
n = 5
# The "house graph"
V = set(range(1, n + 1))
E = {frozenset({1, 2}), frozenset({2, 3}), frozenset({3, 4}), frozenset({4, 1}),
     frozenset({1, 5}), frozenset({2, 5})}

# Define two partitions in P(G, n)
sigma1 = [frozenset({1, 2, 5}), frozenset({3, 4})]
sigma2 = [frozenset({1, 4}), frozenset({2, 3, 5})]

# Compute join and meet
p_join = compute_join(sigma1, sigma2, n)
p_meet = compute_meet(sigma1, sigma2, E)

# Compute ranks
r1 = n - len(sigma1)
r2 = n - len(sigma2)
r_join = n - len(p_join)
r_meet = n - len(p_meet)

print(f"Let G be the house graph with n={n} vertices.")
print(f"Let sigma1 = {format_partition(sigma1)}")
print(f"Let sigma2 = {format_partition(sigma2)}")
print(f"rank(sigma1) = {n} - {len(sigma1)} = {r1}")
print(f"rank(sigma2) = {n} - {len(sigma2)} = {r2}")
print("")
print(f"The join sigma1 V sigma2 = {format_partition(p_join)}")
print(f"rank(join) = {n} - {len(p_join)} = {r_join}")
print("")
print(f"The meet sigma1 ^ sigma2 = {format_partition(p_meet)}")
print(f"rank(meet) = {n} - {len(p_meet)} = {r_meet}")
print("")
print("Checking the semimodular inequality: rank(sigma1) + rank(sigma2) >= rank(join) + rank(meet)")
print(f"{r1} + {r2} >= {r_join} + {r_meet}")
print(f"{r1 + r2} >= {r_join + r_meet}")
