def rotate_90_clockwise(matrix):
    return [list(reversed(col)) for col in zip(*matrix)]

# Test the hypothesis with the provided test input
input_grid = [
    [1, 3, 5],
    [3, 5, 3],
    [8, 8, 3]
]

output_grid = rotate_90_clockwise(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))