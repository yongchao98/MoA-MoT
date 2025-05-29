def is_valid_pattern(grid, r, c):
    rows, cols = len(grid), len(grid[0])
    # Check if position is center of a 3x3 pattern
    if r < 1 or r >= rows-1 or c < 1 or c >= cols-1:
        return False
    if grid[r][c] != 9:
        return False
    # Check surrounding positions
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            if grid[r+dr][c+dc] != 6:
                return False
    return True

def find_connected_group(grid, r, c, visited):
    rows, cols = len(grid), len(grid[0])
    if (r, c) in visited or grid[r][c] not in [6, 9]:
        return set()
    
    group = {(r, c)}
    visited.add((r, c))
    stack = [(r, c)]
    
    while stack:
        curr_r, curr_c = stack.pop()
        # Check all 8 directions
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                new_r, new_c = curr_r + dr, curr_c + dc
                if (0 <= new_r < rows and 0 <= new_c < cols and 
                    (new_r, new_c) not in visited and 
                    grid[new_r][new_c] in [6, 9]):
                    group.add((new_r, new_c))
                    visited.add((new_r, new_c))
                    stack.append((new_r, new_c))
    return group

def transform_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = [[8 for _ in range(cols)] for _ in range(rows)]
    visited = set()
    
    # First find all valid 3x3 patterns
    valid_patterns = set()
    for i in range(rows):
        for j in range(cols):
            if is_valid_pattern(grid, i, j):
                valid_patterns.add((i, j))
    
    # Then find groups connected to bottom
    for j in range(cols):
        if grid[rows-1][j] in [6, 9]:
            group = find_connected_group(grid, rows-1, j, visited)
            # Check if group contains any valid pattern
            has_valid_pattern = any((r, c) in valid_patterns for r, c in group)
            if has_valid_pattern:
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