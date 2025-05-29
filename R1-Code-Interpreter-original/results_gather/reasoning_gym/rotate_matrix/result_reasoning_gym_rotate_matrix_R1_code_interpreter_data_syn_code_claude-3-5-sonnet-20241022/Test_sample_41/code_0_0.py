def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with dimensions reversed
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the new matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Input matrix
matrix = [
    [1, 8, 7, 9, 2, 7],
    [4, 9, 1, 7, 6, 8],
    [8, 3, 7, 1, 0, 4],
    [6, 0, 6, 5, 8, 6],
    [7, 8, 6, 2, 0, 4],
    [9, 9, 3, 1, 7, 5]
]

# Rotate 360 degrees (4 times 90 degrees)
result = matrix
for _ in range(4):
    result = rotate_90_clockwise(result)

# Print the result in the required format
print("Matrix after 360-degree rotation:")
for row in result:
    print(' '.join(map(str, row)))