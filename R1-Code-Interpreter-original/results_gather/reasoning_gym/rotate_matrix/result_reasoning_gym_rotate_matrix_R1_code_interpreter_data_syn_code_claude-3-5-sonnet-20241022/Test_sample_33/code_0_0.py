def rotate_matrix_90(matrix):
    n = len(matrix)
    # Create a new matrix with reversed columns
    return [[matrix[n-1-j][i] for j in range(n)] for i in range(n)]

# Input matrix
matrix = [
    [4, 4, 8, 3, 7, 7, 8, 7],
    [8, 5, 6, 3, 1, 8, 7, 4],
    [3, 7, 9, 0, 4, 1, 2, 8],
    [2, 1, 2, 5, 5, 9, 5, 9],
    [3, 3, 6, 7, 7, 0, 2, 0],
    [5, 4, 5, 4, 0, 2, 6, 1],
    [4, 8, 0, 2, 2, 9, 8, 6],
    [6, 2, 8, 6, 5, 3, 3, 4]
]

# Rotate 720 degrees (4 times 90 degrees)
result = matrix
for _ in range(8):  # 720/90 = 8 rotations
    result = rotate_matrix_90(result)

# Format output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())