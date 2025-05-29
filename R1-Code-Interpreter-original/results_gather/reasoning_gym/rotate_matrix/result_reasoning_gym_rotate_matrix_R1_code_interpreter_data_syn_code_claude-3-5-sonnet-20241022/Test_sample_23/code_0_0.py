def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with reversed rows
    rotated = [[matrix[n-1-j][i] for j in range(n)] for i in range(n)]
    return rotated

# Input matrix
matrix = [
    [3, 2, 4],
    [1, 2, 4],
    [1, 5, 7]
]

# Rotate 8 times (720 degrees)
result = matrix
for _ in range(8):
    result = rotate_90_clockwise(result)

# Format output
output = "\n".join(" ".join(str(num) for num in row) for row in result)
print(f"<<<{output}>>>")