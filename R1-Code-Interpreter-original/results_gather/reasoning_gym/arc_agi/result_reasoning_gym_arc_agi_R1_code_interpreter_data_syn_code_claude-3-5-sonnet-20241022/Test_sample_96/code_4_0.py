def flood_fill_from_edges(grid):
    rows, cols = len(grid), len(grid[0])
    # Keep track of 8's that should remain 8
    keep_8s = set()
    
    def fill(i, j, visited):
        if (i < 0 or i >= rows or j < 0 or j >= cols or
            (i, j) in visited or grid[i][j] != 8):
            return
        
        visited.add((i, j))
        keep_8s.add((i, j))
        
        # Spread to adjacent cells (only up, down, left, right)
        for ni, nj in [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
            fill(ni, nj, visited)
    
    visited = set()
    
    # Start from all edges
    # Top and bottom edges
    for j in range(cols):
        if grid[0][j] == 8:
            fill(0, j, visited)
        if grid[rows-1][j] == 8:
            fill(rows-1, j, visited)
    
    # Left and right edges
    for i in range(rows):
        if grid[i][0] == 8:
            fill(i, 0, visited)
        if grid[i][cols-1] == 8:
            fill(i, cols-1, visited)
    
    return keep_8s

def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = []
    for line in input_grid.strip().split('\n'):
        grid.append([int(x) for x in line.split()])
    
    # Find all 8's that should remain 8
    keep_8s = flood_fill_from_edges(grid)
    
    # Create output grid
    output = []
    for i in range(len(grid)):
        row = []
        for j in range(len(grid[0])):
            if grid[i][j] == 8 and (i,j) not in keep_8s:
                row.append(1)
            else:
                row.append(grid[i][j])
        output.append(row)
    
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