def compare_all_examples(examples):
    results = []
    for i, example in enumerate(examples, 1):
        rows = example.split('\n')
        grid = [list(map(int, row.split())) for row in rows]
        
        # Count clusters
        visited = set()
        clusters = []
        
        def find_cluster(i, j):
            if (i < 0 or i >= len(grid) or j < 0 or j >= len(grid[0]) or 
                (i,j) in visited or grid[i][j] != 2):
                return []
            visited.add((i,j))
            cluster = [(i,j)]
            for ni, nj in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]:
                cluster.extend(find_cluster(ni, nj))
            return cluster
        
        for i in range(len(grid)):
            for j in range(len(grid[0])):
                if grid[i][j] == 2 and (i,j) not in visited:
                    cluster = find_cluster(i, j)
                    if cluster:
                        clusters.append(len(cluster))
        
        results.append({
            'example': i,
            'num_clusters': len(clusters),
            'cluster_sizes': sorted(clusters),
            'total_twos': sum(1 for row in grid for val in row if val == 2)
        })
    return results

# First two examples for comparison
examples = [
"""0 0 3 0 3 3 3
0 0 3 3 3 3 3
3 3 3 3 3 2 2
3 3 0 0 3 2 2
3 3 0 0 3 3 3""",

"""3 3 3 0 0
0 0 3 0 0
0 0 3 3 3
3 3 2 2 3
2 3 2 2 3"""]

comparison = compare_all_examples(examples)
print("Examples comparison:", comparison)