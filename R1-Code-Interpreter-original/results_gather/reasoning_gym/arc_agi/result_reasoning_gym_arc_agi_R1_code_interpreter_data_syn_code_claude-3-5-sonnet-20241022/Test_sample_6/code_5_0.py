def find_groups_from_bottom(grid):
    rows, cols = len(grid), len(grid[0])
    keep = set()  # positions to keep in output
    visited = set()
    
    def dfs(r, c):
        if (r, c) in visited or grid[r][c] not in [6, 9]:
            return
        
        visited.add((r, c))
        keep.add((r, c))
        
        # Check all 8 adjacent positions
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                new_r, new_c = r + dr, c + dc
                if (0 <= new_r < rows and 0 <= new_c < cols and 
                    grid[new_r][new_c] in [6, 9]):
                    dfs(new_r, new_c)
    
    # Start DFS from bottom row
    for j in range(cols):
        if grid[rows-1][j] in [6, 9]:
            dfs(rows-1, j)
    
    return keep

def transform_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = [[8 for _ in range(cols)] for _ in range(rows)]
    
    # Find all positions to keep
    keep_positions = find_groups_from_bottom(grid)
    
    # Build output grid
    for r, c in keep_positions:
        output[r][c] = grid[r][c]
    
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

output = transform_grid(test_input)
print("Final output:")
for row in output:
    print(' '.join(map(str, row)))