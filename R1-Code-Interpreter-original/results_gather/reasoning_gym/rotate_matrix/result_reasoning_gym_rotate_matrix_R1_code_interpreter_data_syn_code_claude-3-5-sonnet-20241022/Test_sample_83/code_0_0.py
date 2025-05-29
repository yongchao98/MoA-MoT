# Input matrix
matrix = '''3 1 8 1 5 8
4 4 5 7 0 8
8 2 2 0 4 2
7 9 3 2 2 2
8 3 2 6 0 3
8 5 9 3 4 9'''

# Convert input string to 2D list
matrix = [list(map(int, row.split())) for row in matrix.split('\n')]

# For 270° clockwise (or 90° counterclockwise):
# 1. Read columns from right to left
# 2. Each column becomes a row
n = len(matrix)
rotated = [[matrix[j][n-1-i] for j in range(n)] for i in range(n)]

# Format output
result = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(result)