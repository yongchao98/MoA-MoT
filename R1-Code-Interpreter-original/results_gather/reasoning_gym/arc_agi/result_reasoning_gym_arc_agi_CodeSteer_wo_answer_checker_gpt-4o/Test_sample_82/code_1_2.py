def find_significant_numbers(grid):
    import numpy as np

    # Convert grid to numpy array for easier manipulation
    grid = np.array(grid)
    rows, cols = grid.shape

    # Step 1: Cluster Detection
    def find_clusters():
        visited = np.zeros_like(grid, dtype=bool)
        clusters = set()

        def dfs(r, c, num):
            stack = [(r, c)]
            cluster_size = 0
            while stack:
                x, y = stack.pop()
                if visited[x, y]:
                    continue
                visited[x, y] = True
                cluster_size += 1
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < rows and 0 <= ny < cols and not visited[nx, ny] and grid[nx, ny] == num:
                        stack.append((nx, ny))
            return cluster_size

        for r in range(rows):
            for c in range(cols):
                if not visited[r, c]:
                    num = grid[r, c]
                    cluster_size = dfs(r, c, num)
                    if cluster_size > 1:  # Consider it a cluster if size > 1
                        clusters.add(num)

        return clusters

    # Step 2: Determine Significant Numbers
    significant_numbers = find_clusters()

    # Step 3: Generate Output
    # Sort significant numbers and repeat them to form the output grid
    significant_numbers = sorted(significant_numbers)
    output_grid = [significant_numbers for _ in range(rows - 1)]

    return output_grid

# Test input grid
input_grid = [
    [1, 8, 1, 1, 1, 1, 1, 1, 1, 1],
    [8, 8, 8, 1, 1, 1, 1, 6, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 6, 6, 1],
    [1, 1, 1, 1, 1, 1, 1, 6, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 3, 3, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 9, 9, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 9, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 9, 1, 1, 1]
]

output_grid = find_significant_numbers(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))