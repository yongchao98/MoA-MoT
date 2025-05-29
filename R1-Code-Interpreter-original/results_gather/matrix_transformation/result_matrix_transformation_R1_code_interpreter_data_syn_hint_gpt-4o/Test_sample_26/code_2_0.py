def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix and then reverse each row
    return [list(reversed(col)) for col in zip(*matrix)]

# Original matrix
matrix = [
    ['4', '2', '9'],
    ['9', '4', '7']
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the result in the required format
formatted_result = ','.join([','.join(row) for row in rotated_matrix])
print(f"<<<{formatted_result}>>>")