def is_valid_3x3(grid, r, c):
    rows, cols = len(grid), len(grid[0])
    if r < 1 or r >= rows-1 or c < 1 or c >= cols-1:
        return False
    if grid[r][c] != 9:
        return False
    # Check all surrounding positions should be 6
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            if grid[r+dr][c+dc] != 6:
                return False
    return True

def find_connected_patterns(grid):
    rows, cols = len(grid), len(grid[0])
    keep = set()
    visited = set()
    
    def dfs_from_pattern(r, c):
        if (r, c) in visited:
            return
        visited.add((r, c))
        
        # If current position is center of valid 3x3
        if is_valid_3x3(grid, r, c):
            # Add the 3x3 pattern to keep set
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    keep.add((r+dr, c+dc))
                    # Continue search from each position of the pattern
                    if grid[r+dr][c+dc] in [6, 9]:
                        dfs_from_pattern(r+dr, c+dc)
        # If current position is a 6 connected to a pattern
        elif grid[r][c] == 6:
            keep.add((r, c))
            # Check neighbors
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    new_r, new_c = r + dr, c + dc
                    if (0 <= new_r < rows and 0 <= new_c < cols and 
                        grid[new_r][new_c] in [6, 9]):
                        dfs_from_pattern(new_r, new_c)
    
    # Start from bottom row
    for j in range(cols):
        if grid[rows-1][j] == 9:
            dfs_from_pattern(rows-1, j)
        elif grid[rows-1][j] == 6:
            # Check if this 6 is part of a valid pattern
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    r, c = rows-1 + dr, j + dc
                    if (0 <= r < rows and 0 <= c < cols and 
                        is_valid_3x3(grid, r, c)):
                        dfs_from_pattern(rows-1, j)
                        break
    
    return keep

def transform_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = [[8 for _ in range(cols)] for _ in range(rows)]
    
    keep_positions = find_connected_patterns(grid)
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