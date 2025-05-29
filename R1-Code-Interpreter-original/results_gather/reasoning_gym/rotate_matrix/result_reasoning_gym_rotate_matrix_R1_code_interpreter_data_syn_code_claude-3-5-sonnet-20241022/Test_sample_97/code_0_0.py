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
    [4, 0, 6, 3, 5, 3, 1, 8, 2, 2],
    [4, 9, 3, 2, 7, 9, 7, 2, 1, 4],
    [1, 4, 1, 6, 0, 3, 9, 2, 1, 3],
    [7, 0, 4, 9, 4, 5, 2, 5, 7, 8],
    [9, 1, 8, 2, 3, 2, 8, 7, 2, 9],
    [5, 4, 9, 5, 9, 7, 5, 5, 2, 8],
    [2, 6, 3, 0, 3, 2, 3, 4, 3, 8],
    [8, 1, 9, 6, 2, 3, 0, 0, 5, 6],
    [2, 0, 2, 8, 2, 9, 4, 6, 4, 6],
    [6, 6, 4, 9, 3, 8, 3, 7, 0, 9]
]

# 450 degrees = 360 + 90 degrees = 1 rotation of 90 degrees
result = rotate_90_clockwise(matrix)

# Format and print the result
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())