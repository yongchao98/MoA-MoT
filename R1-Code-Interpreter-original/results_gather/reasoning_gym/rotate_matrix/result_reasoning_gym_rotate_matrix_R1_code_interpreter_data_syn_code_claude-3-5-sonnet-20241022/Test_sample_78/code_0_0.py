def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with dimensions swapped
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the new matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Input matrix
matrix = [
    [9, 3, 6, 7, 4, 8, 1, 8],
    [6, 8, 5, 0, 1, 8, 1, 2],
    [7, 0, 8, 6, 9, 3, 4, 4],
    [5, 5, 8, 1, 6, 8, 7, 4],
    [2, 3, 4, 6, 7, 3, 2, 7],
    [3, 1, 7, 5, 5, 3, 2, 4],
    [2, 3, 4, 2, 0, 9, 4, 2],
    [0, 7, 6, 3, 3, 7, 8, 7]
]

# Rotate once (90 degrees)
rotated = rotate_90_clockwise(matrix)

# Format and print the output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())