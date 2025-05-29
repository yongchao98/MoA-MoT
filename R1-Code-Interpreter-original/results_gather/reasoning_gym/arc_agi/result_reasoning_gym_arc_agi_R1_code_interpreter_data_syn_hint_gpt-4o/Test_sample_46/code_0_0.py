def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[7] * cols for _ in range(rows)]
    
    def is_contiguous_block(r, c, num):
        # Check if the current position is part of a contiguous block of 'num'
        if r < 0 or r >= rows or c < 0 or c >= cols or input_grid[r][c] != num:
            return False
        # Check surrounding cells
        if (r > 0 and input_grid[r-1][c] == num) or (r < rows-1 and input_grid[r+1][c] == num) or \
           (c > 0 and input_grid[r][c-1] == num) or (c < cols-1 and input_grid[r][c+1] == num):
            return True
        return False
    
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != 7 and is_contiguous_block(r, c, input_grid[r][c]):
                output_grid[r][c] = input_grid[r][c]
    
    return output_grid

# Test input grid
input_grid = [
    [6, 7, 7, 7, 6, 7, 7, 7, 7, 6],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 6, 7, 6, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [6, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 6, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 6, 7],
    [7, 7, 6, 7, 7, 6, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))