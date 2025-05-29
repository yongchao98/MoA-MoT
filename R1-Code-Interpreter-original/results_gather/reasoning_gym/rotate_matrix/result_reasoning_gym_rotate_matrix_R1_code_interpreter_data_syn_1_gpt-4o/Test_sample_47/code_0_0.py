def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Given matrix
matrix = [
    [7, 5, 9, 2],
    [3, 3, 7, 2],
    [1, 5, 0, 3],
    [1, 7, 4, 2]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)