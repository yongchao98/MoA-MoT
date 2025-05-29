def find_connected_cells(grid, start_row, start_col, visited):
    rows, cols = len(grid), len(grid[0])
    connected = set()
    stack = [(start_row, start_col)]
    
    while stack:
        r, c = stack.pop()
        if (r, c) in visited:
            continue
            
        visited.add((r, c))
        connected.add((r, c))
        
        # Check all adjacent cells (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                new_r, new_c = r + dr, c + dc
                if (0 <= new_r < rows and 0 <= new_c < cols and 
                    grid[new_r][new_c] == '7' and 
                    (new_r, new_c) not in visited):
                    stack.append((new_r, new_c))
    
    return connected

def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [row.split() for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Create output grid
    output = [row[:] for row in grid]
    visited = set()
    
    # Find regions that should be converted to 3
    for r in range(rows):
        for c in range(cols):
            # If we find a 7 that's connected to the right edge where there's a 3
            if (grid[r][c] == '7' and c < cols-1 and 
                grid[r][cols-1] == '3' and 
                (r, c) not in visited):
                connected = find_connected_cells(grid, r, c, visited)
                # Convert connected region to 3s
                for cr, cc in connected:
                    output[cr][cc] = '3'
    
    # Format output
    return '\n'.join(' '.join(row) for row in output)

# Test input
test_input = """7 7 7 1 7 7 7 1 7 7 7 1 7
7 7 7 7 7 7 7 7 1 1 1 7 7
7 7 7 1 1 7 1 7 7 7 7 1 7
1 7 7 7 1 7 1 7 7 1 7 7 3
7 7 7 7 7 7 7 7 7 7 7 1 7
1 1 7 7 1 7 1 1 7 7 7 7 7
7 7 1 1 7 7 7 1 7 7 7 1 7
7 1 7 1 7 7 1 7 7 1 7 7 3
7 7 7 7 7 1 7 1 7 7 7 7 7
7 1 7 7 7 7 7 7 1 1 1 7 3
7 1 7 7 7 1 7 7 1 7 1 7 1
1 1 7 7 7 7 7 7 7 1 7 7 1
1 7 7 1 1 7 7 7 7 7 7 7 7"""

print(transform_grid(test_input))