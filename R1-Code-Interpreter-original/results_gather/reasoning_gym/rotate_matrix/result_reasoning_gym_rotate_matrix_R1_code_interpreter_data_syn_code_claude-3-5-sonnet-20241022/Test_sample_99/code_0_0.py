def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with reversed columns
    return [[matrix[n-1-j][i] for j in range(n)] for i in range(n)]

# Input matrix
matrix = [
    [5, 8, 1, 6, 4, 0, 6],
    [1, 5, 9, 9, 1, 2, 3],
    [8, 7, 7, 3, 3, 8, 0],
    [2, 7, 9, 3, 6, 1, 0],
    [7, 4, 6, 1, 2, 1, 2],
    [1, 9, 5, 1, 8, 9, 6],
    [3, 2, 9, 7, 5, 0, 2]
]

# Rotate 3 times (equivalent to 630 degrees)
result = matrix
for _ in range(3):
    result = rotate_90_clockwise(result)

# Format output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())