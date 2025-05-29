# Input matrix
matrix = '''4 9 2 4 3 7
7 5 0 1 2 2
6 0 3 0 6 2
2 4 9 2 8 9
8 9 0 9 9 1
4 2 0 9 3 4'''

# Convert input string to 2D list
matrix = [list(map(int, row.split())) for row in matrix.split('\n')]

n = len(matrix)
rotated = []

# For 270° clockwise:
# Read from right column to left column
# For each column, read from top to bottom
for j in range(n-1, -1, -1):  # columns from right to left
    new_row = []
    for i in range(n):  # rows from top to bottom
        new_row.append(matrix[i][j])
    rotated.append(new_row)

# Convert to output format
result = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{result}>>>")