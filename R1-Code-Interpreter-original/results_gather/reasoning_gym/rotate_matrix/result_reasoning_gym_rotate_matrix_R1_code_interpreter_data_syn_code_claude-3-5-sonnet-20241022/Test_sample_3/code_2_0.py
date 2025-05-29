# Input matrix
matrix = '''5 2 5 2 4 6 4 9 0 8
4 5 9 1 5 0 8 0 1 7
5 0 5 7 0 0 3 5 1 2
5 2 9 1 7 5 1 4 2 9
1 2 8 7 7 5 6 9 6 0
4 9 7 8 8 0 1 5 8 7
1 1 6 6 3 7 6 0 4 9
3 3 9 2 0 1 1 1 2 0
2 7 9 9 0 7 8 2 8 4
8 7 6 0 6 8 9 4 8 8'''

# Convert string to 2D list
matrix = [row.split() for row in matrix.strip().split('\n')]
matrix = [[int(x) for x in row] for row in matrix]

# Create a new matrix for the rotated result
n = len(matrix)
rotated = [[0 for _ in range(n)] for _ in range(n)]

# For 180-degree rotation, each element (i,j) moves to (n-1-i, n-1-j)
for i in range(n):
    for j in range(n):
        rotated[i][j] = matrix[n-1-i][n-1-j]

# Format output
output = '\n'.join(' '.join(str(x) for x in row) for row in rotated)
print(f"<<<{output}>>>")