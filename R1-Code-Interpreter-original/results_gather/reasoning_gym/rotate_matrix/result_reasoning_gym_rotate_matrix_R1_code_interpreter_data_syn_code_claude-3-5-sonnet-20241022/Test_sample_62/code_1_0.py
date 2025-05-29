def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix to store the rotated result
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # For 90-degree clockwise rotation:
    # first row becomes last column
    # second row becomes second-to-last column
    # and so on...
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Input matrix
matrix = [
    [2, 1, 4, 4, 6, 0, 7, 5, 4, 6],
    [9, 4, 8, 0, 6, 0, 6, 5, 7, 5],
    [3, 0, 6, 7, 2, 4, 2, 1, 4, 7],
    [4, 1, 5, 6, 5, 0, 8, 0, 4, 3],
    [0, 7, 5, 3, 4, 2, 0, 2, 7, 9],
    [3, 7, 0, 3, 8, 8, 1, 9, 9, 4],
    [0, 2, 8, 2, 3, 3, 0, 6, 9, 9],
    [9, 2, 3, 2, 3, 7, 2, 8, 6, 3],
    [5, 5, 6, 2, 8, 9, 2, 9, 5, 6],
    [2, 7, 2, 8, 5, 5, 4, 8, 0, 5]
]

# Rotate once for 90 degrees
result = rotate_90_clockwise(matrix)

# Format and print output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())