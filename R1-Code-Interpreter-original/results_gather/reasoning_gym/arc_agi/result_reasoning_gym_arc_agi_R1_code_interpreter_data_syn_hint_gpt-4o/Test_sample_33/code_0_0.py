def generate_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = []

    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != 0:
                # Create a 3x3 block for each non-zero number
                block = [[input_grid[r][c] if (r+i-1 >= 0 and r+i-1 < rows and c+j-1 >= 0 and c+j-1 < cols) else 0 for j in range(3)] for i in range(3)]
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
print(output_grid)