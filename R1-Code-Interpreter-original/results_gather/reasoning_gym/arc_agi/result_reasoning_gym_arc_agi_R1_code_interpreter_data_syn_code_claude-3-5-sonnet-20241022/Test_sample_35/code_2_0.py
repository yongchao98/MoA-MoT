def find_connected_zeros(grid, i, j, visited):
    rows, cols = len(grid), len(grid[0])
    if (i < 0 or i >= rows or j < 0 or j >= cols or 
        grid[i][j] != '0' or (i, j) in visited):
        return set()
    
    visited.add((i, j))
    connected = {(i, j)}
    
    # Check only orthogonal directions (no diagonals)
    directions = [(-1,0), (0,-1), (0,1), (1,0)]
    for di, dj in directions:
        connected.update(find_connected_zeros(grid, i+di, j+dj, visited))
    
    return connected

def should_transform_group(zeros):
    if len(zeros) < 2:
        return False
    
    # Convert to list for easier processing
    zeros = list(zeros)
    
    # Check if forms a 2x2 square
    for i, j in zeros:
        square = {(i,j), (i,j+1), (i+1,j), (i+1,j+1)}
        if square.issubset(zeros):
            return True
    
    # Check if forms a cross or T-shape (3 or more connected orthogonally)
    if len(zeros) >= 3:
        # Check if all zeros are connected orthogonally
        rows = [z[0] for z in zeros]
        cols = [z[1] for z in zeros]
        if max(rows) - min(rows) <= 1 or max(cols) - min(cols) <= 1:
            return True
    
    return False

def transform_grid(input_grid):
    # Convert input string to 2D grid
    grid = [row.split() for row in input_grid.strip().split('\n')]
    output = [row[:] for row in grid]
    
    # Find and process connected zero groups
    visited = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '0' and (i,j) not in visited:
                zeros = find_connected_zeros(grid, i, j, visited)
                if should_transform_group(zeros):
                    for zi, zj in zeros:
                        output[zi][zj] = '5'
    
    return '\n'.join(' '.join(row) for row in output)

# Test input
test_input = """1 7 1 1 1 1 1 1 7 7 1 1 7 7 7 1 1 1
1 1 1 1 1 7 1 1 7 1 1 1 7 0 0 7 0 7
1 7 1 1 7 1 1 7 7 7 1 1 0 7 7 0 0 1
1 7 1 7 7 7 7 7 7 1 1 1 1 1 1 1 1 1
7 7 7 7 1 1 7 1 7 7 1 7 1 1 1 7 7 1
1 1 1 1 7 1 1 7 0 0 7 0 7 1 7 7 7 7
1 1 1 1 7 7 1 7 0 0 0 0 7 1 1 1 1 1
7 1 7 7 7 1 7 7 1 7 1 1 7 1 7 7 1 7
1 1 1 1 7 1 7 7 7 7 7 1 1 1 7 1 1 1
7 1 1 1 1 7 7 1 1 7 7 1 7 1 7 1 7 7
7 7 7 7 1 1 7 1 1 7 0 7 1 7 1 7 7 1
7 1 7 7 7 1 1 1 7 7 7 0 7 1 1 7 7 1
1 7 7 7 1 1 1 7 1 1 7 7 7 1 7 7 1 7
1 7 0 0 7 0 0 7 1 1 1 1 7 7 7 7 7 7
7 0 7 7 0 7 0 1 1 7 1 7 7 7 7 7 7 1
1 7 0 0 7 0 0 1 7 1 7 7 1 1 7 1 7 1
1 7 7 7 1 7 1 1 1 1 1 1 1 1 1 1 1 7"""

print(transform_grid(test_input))