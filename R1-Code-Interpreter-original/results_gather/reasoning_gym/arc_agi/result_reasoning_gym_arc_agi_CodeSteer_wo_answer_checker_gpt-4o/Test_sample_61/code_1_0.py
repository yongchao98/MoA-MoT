def flood_fill(grid, x, y, visited, cluster):
    # Directions for moving in the grid (up, down, left, right)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    stack = [(x, y)]
    visited[x][y] = True
    cluster.append((x, y))
    
    while stack:
        cx, cy = stack.pop()
        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and not visited[nx][ny] and grid[nx][ny] == grid[x][y]:
                visited[nx][ny] = True
                stack.append((nx, ny))
                cluster.append((nx, ny))

def transform_grid(input_grid):
    # Initialize the output grid with the same dimensions as the input grid
    output_grid = [[0] * len(input_grid[0]) for _ in range(len(input_grid))]
    visited = [[False] * len(input_grid[0]) for _ in range(len(input_grid))]
    
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if not visited[i][j]:
                cluster = []
                flood_fill(input_grid, i, j, visited, cluster)
                
                # Determine transformation based on the cluster
                if input_grid[i][j] == 2:
                    # Transform cluster of 2s to 0s
                    for x, y in cluster:
                        output_grid[x][y] = 0
                elif input_grid[i][j] in [9, 8, 4, 5, 1]:
                    # Preserve clusters of 9, 8, 4, 5, 1
                    for x, y in cluster:
                        output_grid[x][y] = input_grid[i][j]
                else:
                    # Handle other numbers (e.g., 7) based on position
                    for x, y in cluster:
                        if x < 3 or y < 3 or x > len(input_grid) - 4 or y > len(input_grid[0]) - 4:
                            output_grid[x][y] = 9
                        else:
                            output_grid[x][y] = 0
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 7, 6, 7],
    [7, 7, 2, 2, 2, 2, 2, 7, 7, 2, 7, 7, 7],
    [7, 7, 2, 7, 7, 7, 2, 7, 7, 2, 2, 7, 7],
    [7, 2, 2, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2],
    [7, 2, 7, 7, 7, 2, 7, 7, 7, 7, 5, 7, 7],
    [7, 2, 7, 9, 7, 2, 7, 7, 7, 7, 7, 7, 7],
    [7, 2, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7],
    [7, 8, 7, 7, 7, 2, 2, 2, 7, 7, 7, 7, 7],
    [7, 2, 7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 2, 2, 7, 7, 7, 3, 7, 7, 7],
    [7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))