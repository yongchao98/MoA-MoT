def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix of the same size
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Input matrix
matrix = [
    [3, 9, 4, 9, 2],
    [5, 2, 9, 4, 9],
    [0, 1, 2, 7, 8],
    [6, 1, 4, 4, 2],
    [7, 2, 7, 2, 0]
]

# Rotate 360 degrees (4 times 90 degrees)
result = matrix
for _ in range(4):
    result = rotate_90_clockwise(result)

# Print the result in the required format
output = "\n".join(" ".join(str(x) for x in row) for row in result)
print(output)