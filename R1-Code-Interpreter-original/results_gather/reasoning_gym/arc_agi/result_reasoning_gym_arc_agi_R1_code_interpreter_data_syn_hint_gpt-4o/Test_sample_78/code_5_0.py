def generate_output_grid(input_list):
    n = len(input_list)
    output_grid = []

    for i in range(n):
        row = [9] * (2 * n - 1)
        row[i] = input_list[i]
        output_grid.append(row)

    return output_grid

input_grid = [6, 7, 4, 5, 9]
output_grid = generate_output_grid(input_grid)

for row in output_grid:
    print(" ".join(map(str, row)))