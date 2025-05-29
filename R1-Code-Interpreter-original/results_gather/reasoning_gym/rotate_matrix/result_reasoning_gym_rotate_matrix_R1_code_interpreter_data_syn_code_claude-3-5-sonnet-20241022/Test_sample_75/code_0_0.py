# Input matrix as string
matrix_str = """9 9 6 7 4
2 6 8 1 4
4 9 0 4 5
3 2 9 3 2
8 2 3 5 1"""

# Convert string to matrix
matrix = [list(map(int, row.split())) for row in matrix_str.strip().split('\n')]

# Get dimensions
n = len(matrix)

# Create rotated matrix by taking columns in reverse order
rotated = [[matrix[n-1-j][i] for j in range(n)] for i in range(n)]

# Format output
output = '\n'.join(' '.join(str(num) for num in row) for row in rotated)

print(f"<<<{output}>>>")