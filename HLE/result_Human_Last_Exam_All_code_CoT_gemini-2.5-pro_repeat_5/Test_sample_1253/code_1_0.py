import collections

class PathCoverSolver:
    """
    Solves the minimum shortest path cover problem using bipartite matching.
    """
    def __init__(self, vertices, edges, start_node):
        # Use sorted list for consistent node-to-index mapping
        self.nodes = sorted(list(vertices))
        self.node_to_idx = {node: i for i, node in enumerate(self.nodes)}
        self.idx_to_node = {i: node for i, node in enumerate(self.nodes)}
        
        self.num_vertices = len(self.nodes)
        self.start_node_str = start_node
        self.start_node_idx = self.node_to_idx[start_node]
        
        self.adj = collections.defaultdict(list)
        for u, v in edges:
            u_idx, v_idx = self.node_to_idx[u], self.node_to_idx[v]
            self.adj[u_idx].append(v_idx)
            self.adj[v_idx].append(u_idx)

    def _bfs(self):
        """Runs BFS from the start node to find shortest distances."""
        distances = [-1] * self.num_vertices
        distances[self.start_node_idx] = 0
        queue = collections.deque([self.start_node_idx])
        
        head = 0
        while queue:
            u = queue.popleft()
            for v in self.adj[u]:
                if distances[v] == -1:
                    distances[v] = distances[u] + 1
                    queue.append(v)
        return distances

    def _dfs_match(self, u, visited_dfs, match, bipartite_adj):
        """DFS helper to find an augmenting path for matching."""
        for v in bipartite_adj[u]:
            if not visited_dfs[v]:
                visited_dfs[v] = True
                if match[v] < 0 or self._dfs_match(match[v], visited_dfs, match, bipartite_adj):
                    match[v] = u
                    return True
        return False

    def _max_bipartite_matching(self, bipartite_adj):
        """Calculates the size of the maximum matching."""
        # match[i] stores the left-partition partner of right-partition node i
        match = [-1] * self.num_vertices
        match_count = 0
        
        # Iterate through all nodes in the left partition (u)
        for u in range(self.num_vertices):
            # visited_dfs must be reset for each attempt to find an augmenting path
            visited_dfs = [False] * self.num_vertices
            if self._dfs_match(u, visited_dfs, match, bipartite_adj):
                match_count += 1
        return match_count

    def solve(self):
        """Executes the full algorithm and prints the result."""
        # 1. Run BFS to get shortest path distances
        distances = self._bfs()
        
        # 2. Build the bipartite graph from the DAG of shortest paths
        bipartite_adj = collections.defaultdict(list)
        for u_idx in range(self.num_vertices):
            for v_idx in self.adj[u_idx]:
                # An edge u -> v exists in the DAG if dist(v) = dist(u) + 1
                if distances[v_idx] == distances[u_idx] + 1:
                    bipartite_adj[u_idx].append(v_idx)
        
        # 3. Find the maximum matching M
        max_matching_size = self._max_bipartite_matching(bipartite_adj)
        
        # 4. The result is |V| - |M|
        min_path_cover = self.num_vertices - max_matching_size
        
        print(f"The problem is solved by finding the minimum path cover in a special DAG.")
        print(f"The size of the minimum path cover is calculated using the formula: |V| - |M|")
        print(f"Where |V| is the number of vertices and |M| is the size of the maximum matching in the corresponding bipartite graph.\n")
        print(f"For the given graph:")
        print(f"Number of vertices |V| = {self.num_vertices}")
        print(f"Size of maximum matching |M| = {max_matching_size}")
        print(f"Final Equation: {self.num_vertices} - {max_matching_size} = {min_path_cover}")
        print(f"The size of the smallest set of paths is: {min_path_cover}")

# --- Main execution using Example 2 ---
V = {'S', '1', '2', '3', '4', '5', '6', '7', '8'}
E = {('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
     ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')}
S = 'S'

solver = PathCoverSolver(V, E, S)
solver.solve()