# Input matrix as string
matrix_str = """0 0 4 1 0 8 2
3 4 6 5 8 4 8
1 9 6 7 6 5 6
3 6 2 1 8 7 3
8 4 2 7 9 5 8
2 0 9 6 7 2 3
6 5 3 7 4 9 0"""

# Convert string to matrix
matrix = [list(map(int, row.split())) for row in matrix_str.strip().split('\n')]

# Get matrix dimensions
n = len(matrix)

# Create rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Perform 270-degree clockwise rotation (equivalent to 90-degree counterclockwise)
for i in range(n):
    for j in range(n):
        rotated[j][n-1-i] = matrix[i][j]

# Format output
output = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{output}>>>")