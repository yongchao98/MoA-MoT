def expand_grid(input_grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Create a copy of the grid to store the output
    output_grid = [row[:] for row in grid]
    
    # Define the expansion pattern
    expansion_patterns = {
        5: [(0, 1), (0, -1)],  # Horizontal expansion
        7: [(1, 0), (-1, 0)],  # Vertical expansion
        9: [(1, 1), (-1, -1)], # Diagonal expansion
        1: [(1, 0), (0, 1)],   # L-shape expansion
        6: [(0, 1), (1, 0)],   # T-shape expansion
        8: [(1, 1), (1, -1)],  # Cross expansion
    }
    
    # Expand each number according to its pattern
    for r in range(rows):
        for c in range(cols):
            num = grid[r][c]
            if num in expansion_patterns:
                for dr, dc in expansion_patterns[num]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        output_grid[nr][nc] = num
    
    # Convert the output grid back to the required format
    return '\n'.join(' '.join(map(str, row)) for row in output_grid)

# Test input grid
test_input = """
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 7 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 7 7 7 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 9 7 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 6 1 1 1 3 3 3 3 3 3 5 9 3 3 3 3 3 3 3
3 3 3 3 3 7 5 3 3 3 3 3 3 3 3 7 6 6 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 3 6 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 6 6 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
"""

# Get the output grid
output_grid = expand_grid(test_input)
print(output_grid)