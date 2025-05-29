# Input matrix
matrix = '''4 9 2 4 3 7
7 5 0 1 2 2
6 0 3 0 6 2
2 4 9 2 8 9
8 9 0 9 9 1
4 2 0 9 3 4'''

# Convert input string to 2D list
matrix = [list(map(int, row.split())) for row in matrix.split('\n')]

# For 270° clockwise (or 90° counterclockwise):
# 1. Read columns from right to left
# 2. Each column becomes a row
n = len(matrix)
rotated = [[matrix[j][n-1-i] for j in range(n)] for i in range(n)]

# Convert to output format
result = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{result}>>>")