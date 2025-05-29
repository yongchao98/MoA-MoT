def rotate_90_clockwise(matrix):
    # Transpose the matrix and then reverse each row
    return [''.join(row) for row in zip(*matrix[::-1])]

# Original matrix
matrix = [
    "FKI",
    "FQY",
    "CIE"
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the result as required
formatted_result = ','.join(','.join(row) for row in rotated_matrix)

# Print the formatted result
print(f"<<<{formatted_result}>>>")