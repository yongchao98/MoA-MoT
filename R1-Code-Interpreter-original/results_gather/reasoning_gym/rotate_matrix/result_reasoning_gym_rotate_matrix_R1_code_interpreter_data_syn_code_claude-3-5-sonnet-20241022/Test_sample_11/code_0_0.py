def rotate_matrix_180(matrix):
    n = len(matrix)
    # Create a new matrix of the same size
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            # For 180° rotation, elements move to opposite positions
            result[n-1-i][n-1-j] = matrix[i][j]
    
    return result

# Input matrix
matrix = [[4, 2],
          [7, 5]]

# Rotate the matrix
rotated = rotate_matrix_180(matrix)

# Print the result in the required format
result = "\n".join(" ".join(str(num) for num in row) for row in rotated)
print(result)