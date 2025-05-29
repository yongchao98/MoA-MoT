def generate_output_grid(input_grid):
    n = len(input_grid)
    output_size = n + input_grid.count(9)
    output_grid = [[9] * output_size for _ in range(output_size)]
    
    for i in range(n):
        output_grid[i][n - i - 1] = input_grid[i]
    
    return output_grid

input_grid = [6, 7, 4, 5, 9]
output_grid = generate_output_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))