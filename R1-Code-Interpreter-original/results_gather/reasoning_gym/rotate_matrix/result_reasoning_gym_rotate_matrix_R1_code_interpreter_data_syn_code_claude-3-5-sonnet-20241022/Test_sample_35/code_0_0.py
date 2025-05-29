def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with dimensions reversed
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the new matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Initial matrix
matrix = [
    [1, 6, 7],
    [6, 4, 9],
    [6, 7, 9]
]

# Rotate 4 times (360 degrees)
result = matrix
for _ in range(4):
    result = rotate_90_clockwise(result)

# Print in required format
output = "\n".join(" ".join(str(x) for x in row) for row in result)
print(output)