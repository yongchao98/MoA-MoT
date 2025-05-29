def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with same dimensions
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    
    return rotated

# Input matrix
matrix = [
    [2, 8, 3, 5, 4, 5, 0, 7, 5],
    [5, 1, 6, 7, 5, 2, 5, 0, 8],
    [1, 2, 5, 0, 6, 2, 2, 6, 4],
    [4, 5, 7, 5, 9, 5, 2, 9, 4],
    [6, 6, 3, 2, 2, 0, 8, 5, 7],
    [0, 8, 6, 1, 4, 0, 6, 2, 1],
    [7, 5, 1, 6, 3, 1, 9, 5, 9],
    [5, 0, 7, 4, 6, 6, 3, 2, 1],
    [6, 7, 1, 4, 7, 6, 0, 5, 3]
]

# Rotate once (90 degrees)
result = rotate_90_clockwise(matrix)

# Format and print the result
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())