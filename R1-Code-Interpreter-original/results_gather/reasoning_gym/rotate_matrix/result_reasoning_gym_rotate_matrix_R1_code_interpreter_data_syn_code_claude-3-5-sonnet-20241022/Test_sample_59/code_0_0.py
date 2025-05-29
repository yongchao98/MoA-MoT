def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with reversed columns
    return [[matrix[n-1-j][i] for j in range(n)] for i in range(n)]

# Input matrix
matrix = [
    [9, 2, 2],
    [6, 7, 7],
    [7, 7, 4]
]

# Rotate 270 degrees clockwise (same as 3 times 90 degrees)
for _ in range(3):
    matrix = rotate_90_clockwise(matrix)

# Format output
result = "\n".join(" ".join(str(num) for num in row) for row in matrix)
print(result)