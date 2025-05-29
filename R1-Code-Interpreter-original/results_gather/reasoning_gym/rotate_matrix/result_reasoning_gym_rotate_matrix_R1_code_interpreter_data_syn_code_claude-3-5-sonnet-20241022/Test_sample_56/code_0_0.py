def rotate_matrix_90(matrix):
    n = len(matrix)
    # Create a new matrix with reversed rows
    return [[matrix[n-1-j][i] for j in range(n)] for i in range(n)]

# Input matrix
matrix = [
    [3, 5, 2, 0, 8, 0, 5, 7],
    [2, 3, 1, 2, 9, 2, 7, 6],
    [5, 6, 2, 4, 2, 1, 7, 5],
    [0, 2, 3, 7, 3, 3, 2, 0],
    [3, 8, 7, 8, 0, 7, 6, 8],
    [6, 7, 8, 1, 7, 2, 6, 0],
    [2, 2, 6, 4, 0, 2, 1, 8],
    [8, 1, 6, 4, 5, 3, 6, 5]
]

# Rotate 720 degrees (2 complete rotations)
result = matrix
for _ in range(2):  # rotate twice for 720 degrees
    result = rotate_matrix_90(rotate_matrix_90(result))

# Print the result in the required format
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())