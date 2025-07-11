import collections

# Hopcroft-Karp algorithm for maximum bipartite matching
# This is a standard implementation, which is required for the chosen algorithm.
def bfs(graph, match, dist, U_nodes):
    queue = collections.deque()
    for u in U_nodes:
        if match[u] == -1:
            dist[u] = 0
            queue.append(u)
        else:
            dist[u] = float('inf')
    dist[-1] = float('inf')
    while queue:
        u = queue.popleft()
        if dist[u] < dist[-1]:
            for v in graph[u]:
                if dist[match[v]] == float('inf'):
                    dist[match[v]] = dist[u] + 1
                    queue.append(match[v])
    return dist[-1] != float('inf')

def dfs(u, graph, match, dist):
    if u != -1:
        for v in graph[u]:
            if dist[match[v]] == dist[u] + 1:
                if dfs(match[v], graph, match, dist):
                    match[v] = u
                    match[u] = v
                    return True
        dist[u] = float('inf')
        return False
    return True

def hopcroft_karp(graph, U_nodes):
    match = {node: -1 for node in list(graph.keys()) + [v for neighbors in graph.values() for v in neighbors]}
    matching = 0
    dist = {}
    while bfs(graph, match, dist, U_nodes):
        for u in U_nodes:
            if match[u] == -1:
                if dfs(u, graph, match, dist):
                    matching += 1
    return matching

def solve():
    # Example 1:
    # V_nodes = {'S', '1', '2', '3', '4'}
    # Edges = {('S', '1'), ('S', '2'), ('1', '3'), ('2', '3'), ('2', '4')}
    # start_node = 'S'
    # |V|=5, Answer=2. |M| must be 3.
    
    # Example 2:
    V_nodes = {'S', '1', '2', '3', '4', '5', '6', '7', '8'}
    Edges = {('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
             ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')}
    start_node = 'S'
    # |V|=9, Answer=3. |M| must be 6.

    # Adjacency list for the undirected graph
    adj = collections.defaultdict(list)
    for u, v in Edges:
        adj[u].append(v)
        adj[v].append(u)

    # Step 1: BFS to find shortest path distances and create layers
    distances = {node: -1 for node in V_nodes}
    queue = collections.deque([start_node])
    distances[start_node] = 0
    
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if distances[v] == -1:
                distances[v] = distances[u] + 1
                queue.append(v)
    
    # Step 2: Construct the shortest-path DAG
    dag_adj = collections.defaultdict(list)
    for u, v_list in adj.items():
        for v in v_list:
            if distances[v] == distances[u] + 1:
                dag_adj[u].append(v)

    # Step 3: Construct the bipartite graph for matching
    # U nodes are all nodes in V, V nodes are also all nodes in V
    # To avoid clashes, we'll label them u_out and v_in internally.
    U_nodes = list(V_nodes)
    bipartite_graph = collections.defaultdict(list)
    v_in_map = {node: f"{node}_in" for node in V_nodes}
    
    all_bipartite_nodes = set(U_nodes)
    for u_out, neighbors in dag_adj.items():
        for v in neighbors:
            v_in = v_in_map[v]
            bipartite_graph[u_out].append(v_in)
            all_bipartite_nodes.add(v_in)
    
    # Step 4: Find max matching
    max_matching_size = hopcroft_karp(bipartite_graph, U_nodes)
    
    # Step 5: Calculate result
    result = len(V_nodes) - max_matching_size

    print(f"Explanation for provided example:")
    print(f"Total number of vertices |V| is {len(V_nodes)}.")
    print(f"The shortest path DAG was constructed.")
    print(f"The maximum matching |M| found in the derived bipartite graph is {max_matching_size}.")
    print(f"The size of the smallest set of paths is |V| - |M|.")
    print(f"{len(V_nodes)} - {max_matching_size} = {result}")

solve()