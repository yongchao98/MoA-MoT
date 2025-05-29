def find_connected_group(grid, r, c, visited):
    rows, cols = len(grid), len(grid[0])
    if (r, c) in visited or grid[r][c] not in [6, 9]:
        return set()
    
    group = {(r, c)}
    visited.add((r, c))
    
    # Check all 8 directions
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            new_r, new_c = r + dr, c + dc
            if (0 <= new_r < rows and 0 <= new_c < cols and 
                (new_r, new_c) not in visited and 
                grid[new_r][new_c] in [6, 9]):
                group.update(find_connected_group(grid, new_r, new_c, visited))
    
    return group

def transform_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = [[8 for _ in range(cols)] for _ in range(rows)]
    visited = set()
    
    # Find all groups starting from each position
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in visited and grid[i][j] in [6, 9]:
                group = find_connected_group(grid, i, j, visited)
                
                # Check if group touches bottom row
                touches_bottom = any(r == rows-1 for r, c in group)
                
                # If group touches bottom, keep it in output
                if touches_bottom:
                    for r, c in group:
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
print("Output grid:")
for row in output:
    print(' '.join(map(str, row)))