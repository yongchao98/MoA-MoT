def extract_2x2_block(grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, row.split())) for row in grid.strip().split('\n')]
    
    # Find the column index of '1' in the middle of the grid
    middle_col = len(grid[0]) // 2
    
    # Iterate over the grid to find the row with '1' in the middle column
    for i in range(len(grid)):
        if grid[i][middle_col] == 1:
            # Extract the 2x2 block from the rows above the identified row
            # Ensure we have enough rows above to extract the block
            if i >= 2:
                block = [[grid[i-2][0], grid[i-2][1]], [grid[i-1][0], grid[i-1][1]]]
                return block

# Test input grid
input_grid = """
8 8 8 8 8 1 8 8 8 8 8
8 8 8 8 8 1 8 8 8 8 8
8 9 9 9 9 1 8 8 8 9 9
8 8 8 8 8 1 8 8 8 8 8
8 8 8 8 8 1 8 8 8 8 8
8 8 8 8 8 1 8 8 8 8 8
8 4 4 4 4 1 8 8 4 4 4
8 8 8 8 8 1 8 8 8 8 8
8 8 8 8 8 1 8 8 8 8 8
"""

# Extract and print the 2x2 block
output_block = extract_2x2_block(input_grid)
for row in output_block:
    print(' '.join(map(str, row)))