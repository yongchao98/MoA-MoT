# Define the original matrix
matrix = [
    [1, 1, 5, 2, 2, 6, 1],
    [2, 9, 0, 9, 8, 1, 9],
    [5, 6, 1, 1, 3, 3, 7],
    [2, 9, 4, 9, 4, 5, 3],
    [7, 1, 1, 0, 2, 9, 1],
    [1, 1, 4, 0, 3, 9, 3],
    [9, 6, 2, 5, 7, 1, 7],
    [6, 2, 8, 2, 8, 1, 4],
    [5, 9, 1, 8, 0, 4, 8]
]

# Function to mirror the matrix along the counterdiagonal
def mirror_counterdiagonal(matrix):
    n = len(matrix)
    m = len(matrix[0])
    mirrored_matrix = [[0] * n for _ in range(m)]
    for i in range(n):
        for j in range(m):
            mirrored_matrix[j][i] = matrix[n - 1 - i][m - 1 - j]
    return mirrored_matrix

# Perform the operations
# Identity transformation (no change)
transformed_matrix = matrix

# Mirror along the counterdiagonal
final_matrix = mirror_counterdiagonal(transformed_matrix)

# Print the final matrix
for row in final_matrix:
    print(' '.join(map(str, row)))