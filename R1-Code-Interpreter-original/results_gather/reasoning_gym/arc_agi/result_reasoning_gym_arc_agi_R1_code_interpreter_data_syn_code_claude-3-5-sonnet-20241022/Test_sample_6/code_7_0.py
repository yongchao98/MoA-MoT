def find_connected_to_bottom(grid):
    rows, cols = len(grid), len(grid[0])
    keep = [[False] * cols for _ in range(rows)]
    visited = set()
    
    def dfs(r, c):
        if (r, c) in visited or grid[r][c] not in [6, 9]:
            return
        visited.add((r, c))
        keep[r][c] = True
        
        # Check all 8 adjacent positions
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                new_r, new_c = r + dr, c + dc
                if (0 <= new_r < rows and 0 <= new_c < cols and 
                    grid[new_r][new_c] in [6, 9]):
                    dfs(new_r, new_c)
    
    # Start from bottom row
    for j in range(cols):
        if grid[rows-1][j] in [6, 9]:
            dfs(rows-1, j)
    
    return keep

def transform_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = [[8] * cols for _ in range(rows)]
    
    # Find positions connected to bottom
    keep = find_connected_to_bottom(grid)
    
    # Keep original values only for connected positions
    for i in range(rows):
        for j in range(cols):
            if keep[i][j]:
                output[i][j] = grid[i][j]
    
    return output

# Test input
test_input = [
    [8,8,8,8,8,8,8,8,6,6,6,8],
    [8,6,6,6,6,6,6,8,6,6,6,8],
    [8,6,9,6,6,6,6,8,8,8,8,8],
    [8,6,6,6,6,9,6,8,6,6,6,6],
    [8,6,6,6,6,6,6,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,6,9,6],
    [8,8,8,8,6,6,6,8,6,6,6,6],
    [8,8,8,8,6,6,9,8,8,8,8,8],
    [8,8,8,8,6,6,6,8,8,8,8,8],
    [6,6,6,8,9,6,6,8,8,8,8,8],
    [9,6,6,8,8,8,8,8,6,9,6,9],
    [6,6,6,8,8,8,8,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,9,6,6]
]

# Print result
result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))