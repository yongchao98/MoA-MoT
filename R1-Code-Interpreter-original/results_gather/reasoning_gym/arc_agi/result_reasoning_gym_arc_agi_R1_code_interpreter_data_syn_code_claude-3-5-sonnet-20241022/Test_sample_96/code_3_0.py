def is_edge(i, j, rows, cols):
    return i == 0 or i == rows-1 or j == 0 or j == cols-1

def find_path_to_edge(grid, start_i, start_j, rows, cols, visited):
    if (start_i, start_j) in visited:
        return False
    
    visited.add((start_i, start_j))
    
    # If we're at an edge and the value is 8, we found a path
    if is_edge(start_i, start_j, rows, cols) and grid[start_i][start_j] == 8:
        return True
    
    # Check all four adjacent cells
    directions = [(-1,0), (1,0), (0,-1), (0,1)]
    for di, dj in directions:
        ni, nj = start_i + di, start_j + dj
        if (0 <= ni < rows and 0 <= nj < cols and 
            grid[ni][nj] == 8 and 
            find_path_to_edge(grid, ni, nj, rows, cols, visited)):
            return True
    
    return False

def get_connected_group(grid, start_i, start_j, rows, cols):
    group = set()
    stack = [(start_i, start_j)]
    
    while stack:
        i, j = stack.pop()
        if (i, j) in group:
            continue
            
        group.add((i, j))
        
        # Check all four adjacent cells
        for ni, nj in [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
            if (0 <= ni < rows and 0 <= nj < cols and 
                grid[ni][nj] == 8 and 
                (ni, nj) not in group):
                stack.append((ni, nj))
    
    return group

def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = []
    for line in input_grid.strip().split('\n'):
        grid.append([int(x) for x in line.split()])
    
    rows = len(grid)
    cols = len(grid[0])
    
    # Create output grid
    output = [[num for num in row] for row in grid]
    processed = set()
    
    # Process each 8 in the grid
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8 and (i, j) not in processed:
                # Get all connected 8's
                group = get_connected_group(grid, i, j, rows, cols)
                processed.update(group)
                
                # Check if this group has a path to edge
                has_path = False
                for gi, gj in group:
                    if find_path_to_edge(grid, gi, gj, rows, cols, set()):
                        has_path = True
                        break
                
                # If no path to edge, convert all 8's in group to 1's
                if not has_path:
                    for gi, gj in group:
                        output[gi][gj] = 1
    
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