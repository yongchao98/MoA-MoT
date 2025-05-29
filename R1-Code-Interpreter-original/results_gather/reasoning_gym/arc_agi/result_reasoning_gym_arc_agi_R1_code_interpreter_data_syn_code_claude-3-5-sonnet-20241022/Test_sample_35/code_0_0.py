def find_connected_zeros(grid, i, j, visited):
    if (i < 0 or i >= len(grid) or j < 0 or j >= len(grid[0]) or 
        grid[i][j] != '0' or (i, j) in visited):
        return set()
    
    visited.add((i, j))
    connected = {(i, j)}
    
    # Check all 8 directions
    directions = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]
    for di, dj in directions:
        connected.update(find_connected_zeros(grid, i+di, j+dj, visited))
    
    return connected

def should_be_five(grid, i, j, zeros):
    # Check if position should be converted to 5 based on patterns
    if len(zeros) < 3:
        return False
        
    # Check cross pattern
    cross = [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]
    cross_count = sum(1 for pos in cross if pos in zeros)
    
    # Check diagonal pattern
    diag = [(i-1,j-1), (i-1,j+1), (i+1,j-1), (i+1,j+1)]
    diag_count = sum(1 for pos in diag if pos in zeros)
    
    return cross_count >= 2 or diag_count >= 2

def transform_grid(input_grid):
    # Convert input string to 2D grid
    grid = [row.split() for row in input_grid.strip().split('\n')]
    
    # Create output grid
    output = [row[:] for row in grid]
    
    # Find all connected zero groups
    visited = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '0' and (i,j) not in visited:
                zeros = find_connected_zeros(grid, i, j, visited)
                # Check each zero if it should be converted to 5
                for zi, zj in zeros:
                    if should_be_five(grid, zi, zj, zeros):
                        output[zi][zj] = '5'
    
    # Convert back to string format
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