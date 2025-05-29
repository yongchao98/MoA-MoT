# Input matrix as string
matrix_str = """8 8 1 2 6 3 4 0
3 1 9 0 1 2 8 4
6 9 6 5 5 1 5 4
9 2 1 8 1 9 1 4
5 1 4 0 5 6 1 7
7 3 3 2 4 3 0 0
6 0 5 5 7 7 9 8
2 3 7 7 5 9 0 4"""

# Convert string to 2D list
matrix = [list(map(int, row.split())) for row in matrix_str.strip().split('\n')]

# Get dimensions
n = len(matrix)

# Create rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Rotate 180 degrees
for i in range(n):
    for j in range(n):
        rotated[i][j] = matrix[n-1-i][n-1-j]

# Format output
output = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(output)