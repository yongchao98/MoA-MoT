def generate_output_grid(input_grid):
    n = len(input_grid)
    output_grid = [[9] * n for _ in range(n)]
    
    for i in range(n):
        output_grid[n - i - 1][i] = input_grid[i]
    
    return output_grid

input_grid = [6, 7, 4, 5, 9]
output_grid = generate_output_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))