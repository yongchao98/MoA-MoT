# Input matrix as string
matrix_str = """8 9 8 8 9 4 1 9
8 7 2 1 9 5 0 0
0 2 9 2 9 9 2 6
7 2 1 2 7 1 6 5
6 0 4 2 2 3 9 2
3 9 2 7 0 6 4 4
4 1 9 8 0 5 6 2
9 9 4 5 1 1 7 7"""

# Convert string to 2D list
matrix = [list(map(int, row.split())) for row in matrix_str.split('\n')]

# Get matrix dimensions
n = len(matrix)

# Create rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Rotate 180 degrees
for i in range(n):
    for j in range(n):
        rotated[i][j] = matrix[n-1-i][n-1-j]

# Format output
output = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{output}>>>")