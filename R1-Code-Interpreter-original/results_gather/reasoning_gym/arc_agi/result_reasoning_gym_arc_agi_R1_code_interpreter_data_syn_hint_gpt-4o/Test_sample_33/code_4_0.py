def generate_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = []

    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != 0:
                # Create a 3x3 block for each non-zero number
                block = []
                for i in range(3):
                    row = []
                    for j in range(3):
                        if 0 <= r+i-1 < rows and 0 <= c+j-1 < cols:
                            row.append(input_grid[r][c])
                        else:
                            row.append(0)
                    block.append(row)
                output_grid.extend(block)
    
    return output_grid

input_grid = [
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 5, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = generate_output_grid(input_grid)
for row in output_grid:
    print(row)