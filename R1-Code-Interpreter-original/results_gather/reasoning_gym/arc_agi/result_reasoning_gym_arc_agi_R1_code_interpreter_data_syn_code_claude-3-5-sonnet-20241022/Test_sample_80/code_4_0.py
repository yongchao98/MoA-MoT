def analyze_two_patterns(grid_str):
    rows = grid_str.split('\n')
    grid = [list(map(int, row.split())) for row in rows]
    
    # Find clusters of 2's (adjacent 2's)
    clusters = []
    visited = set()
    
    def find_cluster(i, j):
        if (i < 0 or i >= len(grid) or j < 0 or j >= len(grid[0]) or 
            (i,j) in visited or grid[i][j] != 2):
            return []
        
        visited.add((i,j))
        cluster = [(i,j)]
        # Check all adjacent positions
        for ni, nj in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]:
            cluster.extend(find_cluster(ni, nj))
        return cluster
    
    # Find all clusters
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == 2 and (i,j) not in visited:
                cluster = find_cluster(i, j)
                if cluster:
                    clusters.append(cluster)
    
    return {
        'num_clusters': len(clusters),
        'cluster_sizes': [len(c) for c in clusters],
        'total_twos': sum(row.count(2) for row in grid)
    }

test_input = """0 0 3 2 3 3 3 3
0 0 3 3 3 3 3 3
3 3 3 3 3 2 2 3
3 2 2 3 3 2 2 3
3 2 2 3 0 3 3 3"""

pattern = analyze_two_patterns(test_input)
print("Test input pattern:", pattern)