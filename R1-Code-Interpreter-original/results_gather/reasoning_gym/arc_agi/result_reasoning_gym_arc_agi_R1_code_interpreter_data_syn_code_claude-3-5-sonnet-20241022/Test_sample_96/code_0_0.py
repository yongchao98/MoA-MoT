def find_connected_8s(grid, i, j, visited):
    if (i < 0 or i >= len(grid) or 
        j < 0 or j >= len(grid[0]) or 
        grid[i][j] != 8 or 
        (i,j) in visited):
        return set()
    
    visited.add((i,j))
    connected = {(i,j)}
    
    # Check all 8 directions
    directions = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]
    for di, dj in directions:
        connected.update(find_connected_8s(grid, i+di, j+dj, visited))
    
    return connected

def reaches_edge(positions, rows, cols):
    for i, j in positions:
        if i == 0 or i == rows-1 or j == 0 or j == cols-1:
            return True
    return False

def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    rows = len(grid)
    cols = len(grid[0])
    
    # Create output grid as a copy
    output = [[num for num in row] for row in grid]
    
    # Find all connected components of 8's
    visited = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8 and (i,j) not in visited:
                connected = find_connected_8s(grid, i, j, visited)
                # If component doesn't reach edge, convert all 8's to 1's
                if not reaches_edge(connected, rows, cols):
                    for pos_i, pos_j in connected:
                        output[pos_i][pos_j] = 1

    # Format output
    result = '\n'.join(' '.join(str(x) for x in row) for row in output)
    print(result)

# Test input
test_input = """8 8 8 0 8 8 4 8 8 8 9 8 2 8 8 5
8 6 5 8 0 8 4 4 8 0 9 7 8 8 8 8
7 8 8 7 3 9 1 8 4 5 7 7 8 8 8 8
0 1 8 8 8 8 8 0 8 8 8 7 4 9 3 8
8 7 8 6 9 5 0 5 8 8 0 4 8 2 9 8
3 8 8 8 4 1 8 8 8 3 8 9 0 1 6 9
2 8 2 8 8 8 8 8 0 8 8 8 8 8 0 8
3 3 1 0 8 8 8 4 8 8 3 3 3 1 8 7
8 8 8 4 8 8 8 8 8 6 8 6 3 4 8 8
2 1 8 3 2 8 8 8 8 8 8 2 8 5 3 8
7 2 8 8 8 8 9 4 1 8 8 3 8 8 8 5
8 2 4 1 7 7 8 2 9 0 5 3 7 8 4 8
8 8 3 8 5 8 9 8 1 3 8 4 8 2 8 0
6 8 8 8 2 8 8 4 8 8 0 8 0 8 9 8
0 7 8 8 8 7 8 8 6 7 1 8 5 8 2 8
1 9 1 5 9 8 8 8 1 7 8 4 8 8 9 8"""

transform_grid(test_input)