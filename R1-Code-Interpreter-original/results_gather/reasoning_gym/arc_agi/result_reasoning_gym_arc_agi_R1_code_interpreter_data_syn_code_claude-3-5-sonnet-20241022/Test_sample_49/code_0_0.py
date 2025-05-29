def find_special_number(grid):
    # Find any number that's not 2 or 5
    special = None
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] not in [2, 5]:
                special = grid[i][j]
                return special, i, j
    return None, -1, -1

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Create a copy
    
    # Find special number and its position
    special_num, special_row, special_col = find_special_number(input_grid)
    
    if special_num is None:
        return output_grid
        
    # Get the block row (special_row // 2 * 2)
    block_row = (special_row // 2) * 2
    
    # Propagate horizontally in the block
    for j in range(cols):
        if input_grid[block_row][j] != 2:
            output_grid[block_row][j] = special_num
            output_grid[block_row + 1][j] = special_num
            
    # Propagate vertically in the same column pattern
    for i in range(0, rows, 3):
        if i + 1 >= rows:
            break
        if input_grid[i][special_col] != 2:
            output_grid[i][special_col] = special_num
            output_grid[i + 1][special_col] = special_num
            if special_col + 1 < cols and input_grid[i][special_col + 1] != 2:
                output_grid[i][special_col + 1] = special_num
                output_grid[i + 1][special_col + 1] = special_num

    return output_grid

# Parse input grid
input_str = """5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
5 2 5 5 2 6 6 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 6 6 2 5 5 2 5 5 2 5 5 2 5 5 2
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2
5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2 5 5 2"""

# Convert input string to grid
input_grid = [[int(x) for x in line.split()] for line in input_str.split('\n')]

# Transform grid
output_grid = transform_grid(input_grid)

# Format output
output_str = '\n'.join(' '.join(str(x) for x in row) for row in output_grid)
print(output_str)