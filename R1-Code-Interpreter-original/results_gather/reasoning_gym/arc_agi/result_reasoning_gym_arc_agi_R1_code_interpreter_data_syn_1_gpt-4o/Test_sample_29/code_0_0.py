def transform_grid(input_grid):
    # Define the size of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Initialize the output grid with the same dimensions
    output_grid = [[0] * cols for _ in range(rows)]
    
    # Function to find clusters using DFS
    def find_cluster(r, c, visited, cluster):
        if (r < 0 or r >= rows or c < 0 or c >= cols or
            visited[r][c] or input_grid[r][c] == 8):
            return
        visited[r][c] = True
        cluster.append((r, c))
        # Explore neighbors
        find_cluster(r+1, c, visited, cluster)
        find_cluster(r-1, c, visited, cluster)
        find_cluster(r, c+1, visited, cluster)
        find_cluster(r, c-1, visited, cluster)
    
    # Find all clusters
    visited = [[False] * cols for _ in range(rows)]
    clusters = []
    for r in range(rows):
        for c in range(cols):
            if not visited[r][c] and input_grid[r][c] != 8:
                cluster = []
                find_cluster(r, c, visited, cluster)
                clusters.append(cluster)
    
    # Apply transformation rule to each cluster
    for cluster in clusters:
        # Determine the transformation value for the cluster
        # For simplicity, let's use the smallest number in the cluster
        cluster_value = min(input_grid[r][c] for r, c in cluster)
        for r, c in cluster:
            output_grid[r][c] = cluster_value
    
    # Fill the rest of the grid with a default value, e.g., 3
    for r in range(rows):
        for c in range(cols):
            if output_grid[r][c] == 0:
                output_grid[r][c] = 3
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 8, 8],
    [8, 8, 8, 0, 8, 8, 8, 2, 2, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 8],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8],
    [8, 8, 8, 8, 8, 2, 8, 8, 8, 8, 8, 2, 8],
    [8, 8, 8, 8, 8, 2, 2, 2, 8, 8, 8, 6, 8],
    [8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 2, 8],
    [8, 8, 8, 8, 8, 8, 8, 2, 8, 4, 8, 2, 8],
    [8, 8, 9, 8, 8, 8, 8, 2, 8, 8, 8, 2, 8],
    [2, 2, 2, 8, 8, 8, 2, 2, 8, 8, 2, 2, 8],
    [8, 8, 2, 2, 8, 8, 2, 8, 8, 8, 2, 8, 8],
    [8, 8, 8, 2, 8, 8, 2, 2, 2, 2, 2, 8, 8],
    [8, 7, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))