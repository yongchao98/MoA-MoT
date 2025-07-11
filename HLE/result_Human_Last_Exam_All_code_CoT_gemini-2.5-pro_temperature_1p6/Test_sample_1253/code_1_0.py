import collections

def solve_path_cover():
    """
    Solves the minimum shortest path cover problem for a given graph.
    The graph from Example 1 is hardcoded for demonstration.
    """
    # Example 1:
    # V = {'S', '1', '2', '3', '4'}
    # E = {('S', '1'), ('S', '2'), ('1', '3'), ('2', '3'), ('2', '4')}
    # S_node = 'S'
    # Expected output: 2

    # Example 2:
    # V = {'S', '1', '2', '3', '4', '5', '6', '7', '8'}
    # E = {('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
    #      ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')}
    # S_node = 'S'
    # Expected output: 3
    
    # Using Example 1 for demonstration
    V_list = ['S', '1', '2', '3', '4']
    E = [('S', '1'), ('S', '2'), ('1', '3'), ('2', '3'), ('2', '4')]
    S_node = 'S'

    # Step 0: Create adjacency list for the undirected graph
    adj = collections.defaultdict(list)
    for u, v in E:
        adj[u].append(v)
        adj[v].append(u)

    # Step 1: BFS to find shortest distances from S
    dist = {node: -1 for node in V_list}
    q = collections.deque([S_node])
    dist[S_node] = 0
    
    head = 0
    while head < len(q):
        u = q.popleft()
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)

    # Step 2: Build the Shortest-Path DAG
    dag = collections.defaultdict(list)
    for u, v_neighbor in E:
        # Check both directions since original graph is undirected
        if dist[v_neighbor] == dist[u] + 1:
            dag[u].append(v_neighbor)
        elif dist[u] == dist[v_neighbor] + 1:
            dag[v_neighbor].append(u)
    
    # Step 3: Find Maximum Bipartite Matching
    # The left partition is U=V_list, the right partition is V=V_list.
    # An edge from u in U to v in V exists if u->v is in the DAG.
    
    # match_R[v] stores the node u in the left partition matched with v.
    match_R = {node: None for node in V_list}

    def dfs_match(u, visited_dfs, current_dag, current_match_R):
        """DFS to find an augmenting path."""
        for v in current_dag.get(u, []):
            if not visited_dfs[v]:
                visited_dfs[v] = True
                # If v is not matched or its partner can find a new match
                if current_match_R[v] is None or dfs_match(current_match_R[v], visited_dfs, current_dag, current_match_R):
                    current_match_R[v] = u
                    return True
        return False

    matching_size = 0
    # Iterate through all nodes in the left partition
    for u in V_list:
        # For each u, try to find an augmenting path
        # visited_dfs must be reset for each starting node u
        visited_dfs = {node: False for node in V_list}
        if dfs_match(u, visited_dfs, dag, match_R):
            matching_size += 1

    # Step 4: Calculate the result using Konig's theorem
    num_vertices = len(V_list)
    min_path_cover_size = num_vertices - matching_size

    print("The size of the minimum path cover is found by the formula |V| - |M|:")
    print(f"Number of vertices |V| = {num_vertices}")
    print(f"Size of maximum matching |M| = {matching_size}")
    print(f"Result: {num_vertices} - {matching_size} = {min_path_cover_size}")


solve_path_cover()