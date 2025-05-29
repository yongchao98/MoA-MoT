def rotate_matrix(matrix):
    n = len(matrix)
    # Create a new matrix of the same size
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # For clockwise rotation:
    # first row becomes last column
    # second row becomes second-to-last column
    # and so on...
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Initial matrix
matrix = [
    [9, 5, 2, 8],
    [9, 2, 6, 4],
    [9, 6, 9, 9],
    [0, 8, 8, 9]
]

# Rotate 720 degrees (2 complete rotations)
rotations = 2 * 4  # 4 rotations = 360 degrees
current = matrix
for _ in range(rotations):
    current = rotate_matrix(current)

# Print the result in the required format
result = ""
for row in current:
    result += " ".join(map(str, row)) + "\n"
print(result.strip())