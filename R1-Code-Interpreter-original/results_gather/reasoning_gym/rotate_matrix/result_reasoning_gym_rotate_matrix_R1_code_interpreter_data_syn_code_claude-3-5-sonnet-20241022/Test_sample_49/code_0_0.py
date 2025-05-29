def rotate_180(matrix):
    n = len(matrix)
    # Create a new matrix to store the result
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Rotate 180 degrees
    for i in range(n):
        for j in range(n):
            result[n-1-i][n-1-j] = matrix[i][j]
    
    return result

# Input matrix
matrix = [
    [7, 5, 8, 8, 6, 7, 0, 2, 2, 7],
    [9, 0, 5, 4, 5, 3, 4, 2, 9, 1],
    [3, 8, 4, 9, 8, 2, 5, 2, 6, 3],
    [7, 3, 7, 0, 8, 6, 5, 4, 5, 4],
    [9, 5, 1, 2, 0, 6, 4, 4, 2, 3],
    [7, 1, 3, 3, 1, 2, 9, 4, 6, 0],
    [0, 6, 8, 2, 2, 1, 6, 6, 3, 0],
    [9, 5, 0, 3, 5, 7, 6, 4, 1, 6],
    [0, 6, 3, 1, 0, 8, 7, 3, 5, 4],
    [9, 9, 0, 7, 8, 6, 8, 5, 3, 4]
]

# Rotate the matrix
rotated = rotate_180(matrix)

# Print the result in the required format
result = ""
for row in rotated:
    result += " ".join(map(str, row)) + "\n"
print(result.strip())