def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix of same size
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Initial matrix
matrix = [[9, 9],
          [7, 3]]

# Rotate 10 times (900 degrees)
result = matrix
for _ in range(10):
    result = rotate_90_clockwise(result)

# Print the result in required format
print('\n'.join(' '.join(str(x) for x in row) for row in result))