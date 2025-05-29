def get_connected_component(grid, start_i, start_j):
    rows, cols = len(grid), len(grid[0])
    component = set()
    queue = [(start_i, start_j)]
    
    while queue:
        i, j = queue.pop(0)
        if (i, j) in component:
            continue
            
        component.add((i, j))
        
        # Check four adjacent cells
        for ni, nj in [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
            if (0 <= ni < rows and 0 <= nj < cols and 
                grid[ni][nj] == 8 and 
                (ni, nj) not in component):
                queue.append((ni, nj))
    
    return component

def is_trapped(component, rows, cols):
    # Check if any position in the component touches an edge
    for i, j in component:
        if i == 0 or i == rows-1 or j == 0 or j == cols-1:
            return False
    return True

def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = []
    for line in input_grid.strip().split('\n'):
        grid.append([int(x) for x in line.split()])
    
    rows, cols = len(grid), len(grid[0])
    processed = set()
    output = [[x for x in row] for row in grid]
    
    # Process each unprocessed 8
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8 and (i, j) not in processed:
                # Get the connected component containing this 8
                component = get_connected_component(grid, i, j)
                processed.update(component)
                
                # If the component is trapped, convert all its 8's to 1's
                if is_trapped(component, rows, cols):
                    for ci, cj in component:
                        output[ci][cj] = 1
    
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