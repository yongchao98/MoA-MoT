def print_matrix(matrix):
    for row in matrix:
        print(' '.join(map(str, row)))

# Input matrix
matrix = [
    [8, 5, 7, 0, 8],
    [3, 9, 3, 4, 7],
    [3, 5, 7, 1, 4],
    [3, 2, 1, 7, 8],
    [5, 3, 7, 5, 3]
]

# 900° = 720° + 180° = (360° × 2) + 180°
# 900° ÷ 90° = 10 rotations
# After every 4 rotations, we get back to original
# 10 rotations = 4 × 2 + 2 = 2 rotations

def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix
    result = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    return result

# Perform 2 rotations (equivalent to 900 degrees)
for _ in range(2):
    matrix = rotate_90_clockwise(matrix)

# Print result
result = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(result)