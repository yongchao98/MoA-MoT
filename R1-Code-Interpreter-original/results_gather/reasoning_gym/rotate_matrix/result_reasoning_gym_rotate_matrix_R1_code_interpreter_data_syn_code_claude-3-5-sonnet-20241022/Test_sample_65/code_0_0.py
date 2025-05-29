# Input matrix
matrix = [
    [6, 5, 0, 9, 0],
    [8, 0, 0, 4, 5],
    [0, 1, 7, 7, 1],
    [5, 3, 4, 2, 4],
    [6, 8, 3, 2, 1]
]

def rotate_270_clockwise(matrix):
    n = len(matrix)
    # For 270° clockwise rotation, we can transpose and then reverse columns
    # First, create a new matrix with swapped rows and columns (transpose)
    rotated = [[matrix[j][i] for j in range(n)] for i in range(n)]
    # Then reverse each row
    rotated = [row[::-1] for row in rotated]
    return rotated

# Rotate the matrix
result = rotate_270_clockwise(matrix)

# Format output
output = "\n".join(" ".join(str(num) for num in row) for row in result)
print(f"<<<{output}>>>")