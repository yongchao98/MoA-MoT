def rotate_180(matrix):
    # Reverse each row
    for i in range(len(matrix)):
        matrix[i] = matrix[i][::-1]
    
    # Reverse the order of rows
    matrix = matrix[::-1]
    
    return matrix

# Input matrix
matrix = [
    [0, 7, 0, 8, 5, 3],
    [7, 7, 8, 5, 9, 2],
    [8, 4, 8, 7, 6, 5],
    [9, 3, 2, 7, 0, 8],
    [0, 2, 0, 7, 8, 4],
    [8, 1, 0, 3, 3, 6]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)