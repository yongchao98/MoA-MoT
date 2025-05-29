def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed = list(zip(*matrix))
    # Reverse each row in the transposed matrix
    rotated = [list(row)[::-1] for row in transposed]
    return rotated

# Original matrix
matrix = [
    ['9', 'H', 'B', 'X'],
    ['S', 'X', 'T', '6'],
    ['U', 'V', 'N', '6']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output as requested
formatted_output = ','.join([','.join(row) for row in rotated_matrix])
print(f"<<<{formatted_output}>>>")