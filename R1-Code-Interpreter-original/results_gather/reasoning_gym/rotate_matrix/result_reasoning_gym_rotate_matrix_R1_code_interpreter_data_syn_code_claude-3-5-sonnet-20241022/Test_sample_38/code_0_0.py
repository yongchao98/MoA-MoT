def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with dimensions swapped
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    
    return rotated

# Input matrix
matrix = [
    [4, 7, 7],
    [2, 5, 1],
    [2, 9, 3]
]

# Rotate 10 times (900 degrees)
for _ in range(10):
    matrix = rotate_90_clockwise(matrix)

# Format output
output = "\n".join(" ".join(str(num) for num in row) for row in matrix)
print(output)