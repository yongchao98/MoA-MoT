import collections

def solve_min_path_cover():
    """
    Solves the minimum shortest path cover problem for a given graph.
    This implementation uses the algorithm for Example 2.
    """
    # Graph definition for Example 2
    # V = {S, 1, 2, 3, 4, 5, 6, 7, 8}, S is mapped to 0
    num_vertices = 9
    edges = [
        (0, 1), (0, 2), (0, 3), (1, 4), (2, 4), 
        (2, 5), (3, 5), (4, 6), (4, 7), (5, 8)
    ]
    source = 0
    
    print("Solving for Example 2:")
    print("V = {S, 1, 2, 3, 4, 5, 6, 7, 8}, S=0")
    print("E = {(S, 1), (S, 2), (S, 3), (1, 4), (2, 4), (2, 5), (3, 5), (4, 6), (4, 7), (5, 8)}")
    
    # Adjacency list for the original undirected graph
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # --- Step 1: Create the Shortest Path DAG using BFS ---
    distances = [-1] * num_vertices
    distances[source] = 0
    queue = collections.deque([source])
    
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if distances[v] == -1:
                distances[v] = distances[u] + 1
                queue.append(v)

    dag_adj = collections.defaultdict(list)
    for u in range(num_vertices):
        for v in adj[u]:
            if distances[u] != -1 and distances[v] == distances[u] + 1:
                dag_adj[u].append(v)
    
    # --- Step 2 & 3: Find Maximum Bipartite Matching on the DAG ---
    
    # match[v] = u means vertex v is matched with vertex u (u -> v)
    match = [-1] * num_vertices
    
    def dfs_for_matching(u, visited, current_match, current_dag):
        """DFS to find an augmenting path."""
        for v in current_dag[u]:
            if not visited[v]:
                visited[v] = True
                # If v is unmatched, or its current match can be re-matched
                if current_match[v] < 0 or dfs_for_matching(current_match[v], visited, current_match, current_dag):
                    current_match[v] = u
                    return True
        return False

    matching_size = 0
    for u in range(num_vertices):
        visited = [False] * num_vertices
        if dfs_for_matching(u, visited, match, dag_adj):
            matching_size += 1

    # --- Step 4: Calculate result using |V| - |M| ---
    min_paths = num_vertices - matching_size
    
    print("\nCalculation based on the algorithm (|V| - |M|):")
    print(f"Total number of vertices |V| = {num_vertices}")
    print(f"Size of the maximum matching |M| = {matching_size}")
    print(f"Size of the smallest set of paths = |V| - |M| = {num_vertices} - {matching_size} = {min_paths}")

if __name__ == '__main__':
    solve_min_path_cover()