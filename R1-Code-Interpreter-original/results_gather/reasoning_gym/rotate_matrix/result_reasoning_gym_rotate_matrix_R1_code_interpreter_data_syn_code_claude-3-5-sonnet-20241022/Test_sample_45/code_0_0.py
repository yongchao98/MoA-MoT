def rotate_matrix_90(matrix):
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
    [5, 3, 3, 3, 1, 4],
    [7, 8, 9, 3, 2, 8],
    [8, 0, 4, 0, 7, 2],
    [3, 8, 6, 2, 5, 5],
    [1, 1, 8, 4, 8, 7],
    [3, 3, 3, 1, 1, 3]
]

# Rotate 720 degrees (8 times 90 degrees)
result = matrix
for _ in range(8):  # 8 rotations of 90 degrees = 720 degrees
    result = rotate_matrix_90(result)

# Print the result in the required format
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())