# Original matrix
matrix = [
    [5, 9],
    [9, 0]
]

# Function to rotate the matrix 90 degrees counterclockwise
def rotate_90_counterclockwise(mat):
    # Transpose the matrix and then reverse each row
    return [list(row) for row in zip(*mat)][::-1]

# Rotate the matrix
rotated_matrix = rotate_90_counterclockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)