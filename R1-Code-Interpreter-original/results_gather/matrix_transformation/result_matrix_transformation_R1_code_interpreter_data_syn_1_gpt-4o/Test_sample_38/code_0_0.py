def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the rotated matrix
    rotated_matrix = [list(reversed(row)) for row in transposed_matrix]
    return rotated_matrix

# Original matrix
matrix = [
    ["4", "1", "3"],
    ["4", "7", "1"]
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output as a list of values and symbols
formatted_output = ",".join([",".join(row) for row in rotated_matrix])
print(formatted_output)