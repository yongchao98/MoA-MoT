def find_connected_groups(grid, i, j, visited, target):
    if (i < 0 or i >= len(grid) or j < 0 or j >= len(grid[0]) or 
        visited[i][j] or grid[i][j] != target or target == 3):
        return []
    
    connected = [(i, j)]
    visited[i][j] = True
    
    # Check all 8 directions
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            if di == 0 and dj == 0:
                continue
            ni, nj = i + di, j + dj
            if (0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and 
                not visited[ni][nj] and grid[ni][nj] == target):
                connected.extend(find_connected_groups(grid, ni, nj, visited, target))
    
    return connected

def process_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # First pass: find initial groups
    visited = [[False] * cols for _ in range(rows)]
    groups_to_convert = []
    
    for i in range(rows):
        for j in range(cols):
            if not visited[i][j] and input_grid[i][j] != 3:
                group = find_connected_groups(input_grid, i, j, visited, input_grid[i][j])
                if len(group) >= 3:
                    groups_to_convert.append(group)
    
    # Convert groups to 5's
    for group in groups_to_convert:
        for x, y in group:
            output_grid[x][y] = 5
    
    # Print the result
    result = []
    for row in output_grid:
        result.append(' '.join(map(str, row)))
    print('\n'.join(result))

# Test input
test_input = [
    [3,9,3,3,0,3,2,4,3,3,6,3,3,3,6],
    [3,3,3,5,3,3,7,3,3,2,6,3,3,0,0],
    [3,6,3,3,3,6,0,3,2,7,3,2,3,3,3],
    [7,3,7,3,0,0,3,5,0,3,3,3,3,3,3],
    [3,4,2,9,3,0,3,3,2,7,3,5,2,3,3],
    [3,7,3,3,5,3,3,3,0,2,4,3,3,3,3],
    [0,3,0,3,3,3,3,4,7,3,3,2,3,0,2],
    [3,3,9,9,3,6,7,3,0,7,3,4,0,4,3],
    [7,7,7,0,3,2,7,3,5,3,9,7,3,5,3],
    [3,3,3,3,9,3,6,3,0,2,3,0,4,3,3],
    [3,7,5,3,3,3,6,3,0,3,3,7,3,4,9],
    [3,3,3,3,3,2,2,7,3,3,4,2,3,3,9],
    [3,3,9,3,2,3,3,3,6,2,2,3,4,3,3],
    [3,3,7,7,0,2,3,2,3,3,6,2,3,5,3],
    [4,2,9,2,6,7,7,3,7,6,2,7,3,4,0]
]

process_grid(test_input)