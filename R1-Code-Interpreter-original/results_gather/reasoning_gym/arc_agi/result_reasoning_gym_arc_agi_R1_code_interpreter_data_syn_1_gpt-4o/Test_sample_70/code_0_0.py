def generate_output_grid(input_grid):
    output_grid = []
    for i in range(3):
        row = input_grid[i] + input_grid[i][::-1]
        output_grid.append(row)
    output_grid.append(output_grid[2])
    output_grid.append(output_grid[1])
    output_grid.append(output_grid[0])
    return output_grid

input_grid = [
    [9, 9, 5],
    [7, 7, 4],
    [7, 7, 7]
]

output_grid = generate_output_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))