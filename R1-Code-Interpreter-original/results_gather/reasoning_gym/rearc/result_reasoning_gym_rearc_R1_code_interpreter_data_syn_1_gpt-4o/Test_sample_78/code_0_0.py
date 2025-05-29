def apply_rule(input_grid):
    # Convert the input grid into a list of lists
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    
    # Check the last element of the last row
    if grid[-1][-1] != 2:
        grid[-1][-1] = 2
    
    # Convert the grid back to the required output format
    output_grid = '\n'.join(' '.join(map(str, row)) for row in grid)
    return output_grid

# Test input grid
input_grid = """
7 2 3 6 3 8 1
1 6 8 2 1 7 0
1 3 2 5 7 2 3
2 3 0 7 8 1 4
9 6 7 4 1 9 7
0 9 2 3 3 2 4
7 2 6 3 6 0 9
"""

# Apply the rule to the test input grid
output_grid = apply_rule(input_grid)
print(output_grid)