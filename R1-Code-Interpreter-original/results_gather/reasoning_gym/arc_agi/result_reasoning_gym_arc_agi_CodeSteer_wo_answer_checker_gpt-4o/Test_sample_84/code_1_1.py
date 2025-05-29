def find_output_length(grid):
    # Check the number of '9's in the main diagonal
    n = len(grid)
    count_nines = 0
    for i in range(n):
        if grid[i][i] == 9:
            count_nines += 1
    return count_nines

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 6, 6, 6, 9, 9, 9],
    [9, 9, 9, 6, 9, 6, 6, 6, 9],
    [9, 6, 6, 6, 9, 9, 9, 6, 6],
    [6, 6, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 6, 6, 6, 9, 9],
    [6, 6, 9, 9, 6, 9, 6, 6, 6],
    [9, 6, 6, 6, 6, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9]
]

output_length = find_output_length(input_grid)
output_grid = [9] * output_length
print("<<<" + " ".join(map(str, output_grid)) + ">>>")