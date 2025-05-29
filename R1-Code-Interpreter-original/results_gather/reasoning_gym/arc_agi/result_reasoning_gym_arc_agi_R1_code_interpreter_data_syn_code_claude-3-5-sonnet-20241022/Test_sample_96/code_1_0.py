def can_reach_edge(grid, start_i, start_j, rows, cols):
    if start_i == 0 or start_i == rows-1 or start_j == 0 or start_j == cols-1:
        return True
        
    visited = set()
    stack = [(start_i, start_j)]
    
    while stack:
        i, j = stack.pop()
        if (i, j) in visited:
            continue
            
        visited.add((i, j))
        
        # Check four adjacent positions
        for ni, nj in [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
            if 0 <= ni < rows and 0 <= nj < cols and grid[ni][nj] == 8:
                # If we reach an edge, this path works
                if ni == 0 or ni == rows-1 or nj == 0 or nj == cols-1:
                    return True
                # If we haven't visited this position yet, add it to stack
                if (ni, nj) not in visited:
                    stack.append((ni, nj))
    
    return False

def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = []
    for line in input_grid.strip().split('\n'):
        grid.append([int(x) for x in line.split()])
    
    rows, cols = len(grid), len(grid[0])
    output = [[x for x in row] for row in grid]
    
    # Check each 8 individually
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8:
                # If this 8 cannot reach an edge through a path of 8's, convert to 1
                if not can_reach_edge(grid, i, j, rows, cols):
                    output[i][j] = 1
    
    # Format output
    for row in output:
        print(' '.join(str(x) for x in row))

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