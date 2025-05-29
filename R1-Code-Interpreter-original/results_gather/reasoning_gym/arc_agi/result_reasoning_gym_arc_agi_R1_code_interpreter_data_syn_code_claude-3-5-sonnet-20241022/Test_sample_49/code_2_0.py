def transform_grid(input_grid):
    rows, cols = len(input_grid), len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # Find special number (not 2 or 5)
    special_num = None
    special_row = special_col = -1
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] not in [2, 5]:
                special_num = input_grid[i][j]
                special_row = i
                special_col = j
                break
        if special_num:
            break
            
    if not special_num:
        return output_grid
        
    # Get block row and column
    block_row = (special_row // 3) * 3
    block_col = (special_col // 3) * 3
    
    # Process horizontal propagation
    for j in range(0, cols, 3):
        if j + 2 >= cols:
            break
        # Replace the 2×2 block in the same row block
        if input_grid[block_row][j] != 2:
            output_grid[block_row][j] = special_num
            output_grid[block_row][j+1] = special_num
            output_grid[block_row+1][j] = special_num
            output_grid[block_row+1][j+1] = special_num
            
    # Process vertical propagation
    for i in range(0, rows, 3):
        if i + 2 >= rows:
            break
        # Replace the 2×2 block in the same column pattern
        if input_grid[i][block_col] != 2:
            output_grid[i][block_col] = special_num
            output_grid[i][block_col+1] = special_num
            output_grid[i+1][block_col] = special_num
            output_grid[i+1][block_col+1] = special_num
            
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

input_grid = [[int(x) for x in line.split()] for line in input_str.split('\n')]
output_grid = transform_grid(input_grid)
output_str = '\n'.join(' '.join(str(x) for x in row) for row in output_grid)
print(output_str)