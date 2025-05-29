# Input matrix
matrix = '''6 3 9 2 3 3 2 1
4 6 8 8 1 9 7 3
6 9 9 3 5 6 3 2
7 2 9 5 1 8 1 7
5 6 0 4 8 0 6 3
5 8 1 7 8 3 9 7
0 5 2 3 2 7 9 2
8 2 0 1 0 2 7 6'''

# Convert string to 2D list
matrix = [list(map(int, row.split())) for row in matrix.strip().split('\n')]

# Get dimensions
n = len(matrix)

# Create new rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Rotate 90 degrees clockwise
for i in range(n):
    for j in range(n):
        rotated[j][n-1-i] = matrix[i][j]

# Format output
output = ''
for row in rotated:
    output += ' '.join(map(str, row)) + '\n'
output = output.strip()

print(f"<<<{output}>>>")