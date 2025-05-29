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
    [5, 7, 3, 6, 6],
    [9, 2, 3, 6, 8],
    [3, 2, 7, 6, 6],
    [4, 6, 3, 5, 7],
    [0, 5, 3, 9, 7]
]

# Rotate the matrix
rotated = rotate_180(matrix)

# Print the result in the required format
result = ""
for row in rotated:
    result += " ".join(map(str, row)) + "\n"
print(result.strip())