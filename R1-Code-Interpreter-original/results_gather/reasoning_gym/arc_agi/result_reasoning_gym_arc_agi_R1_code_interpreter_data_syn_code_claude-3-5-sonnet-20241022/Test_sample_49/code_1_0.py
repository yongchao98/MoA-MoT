def find_special_number(grid):
    # Find any number that's not 2 or 5
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] not in [2, 5]:
                return grid[i][j], i, j
    return None, -1, -1

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # Find special number and its position
    special_num, special_row, special_col = find_special_number(input_grid)
    
    if special_num is None:
        return output_grid

    # Find the column pattern (the 2×2 block position)
    block_col = (special_col // 3) * 3
    
    # For each pair of rows (excluding separator rows)
    for i in range(0, rows, 3):
        if i + 1 >= rows:
            break
        # Replace the corresponding 2×2 block with special number
        if input_grid[i][special_col] != 2:  # Skip separator columns
            output_grid[i][block_col + 1] = special_num
            output_grid[i + 1][block_col + 1] = special_num
            if block_col + 2 < cols:
                output_grid[i][block_col + 2] = special_num
                output_grid[i + 1][block_col + 2] = special_num

    # Process the row block where special number appears
    special_block_row = (special_row // 3) * 3
    for j in range(0, cols, 3):
        if j + 2 >= cols:
            break
        if input_grid[special_row][j] != 2:  # Skip separator columns
            output_grid[special_block_row][j] = special_num
            output_grid[special_block_row][j + 1] = special_num
            output_grid[special_block_row + 1][j] = special_num
            output_grid[special_block_row + 1][j + 1] = special_num

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