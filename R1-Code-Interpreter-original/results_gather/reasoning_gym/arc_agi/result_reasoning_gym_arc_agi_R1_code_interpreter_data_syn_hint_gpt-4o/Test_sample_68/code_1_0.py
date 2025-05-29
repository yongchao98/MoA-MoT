def transform_grid(matrix):
    # Rotate 90 degrees clockwise
    rotated = [list(reversed(col)) for col in zip(*matrix)]
    # Reverse each row
    transformed = [row[::-1] for row in rotated]
    return transformed

# Test the hypothesis with the provided test input
input_grid = [
    [1, 3, 5],
    [3, 5, 3],
    [8, 8, 3]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))