# Input matrix
matrix = '''4 9 2 4 3 7
7 5 0 1 2 2
6 0 3 0 6 2
2 4 9 2 8 9
8 9 0 9 9 1
4 2 0 9 3 4'''

# Convert input string to 2D list
matrix = [list(map(int, row.split())) for row in matrix.split('\n')]

# For 270° clockwise rotation:
# We can do this by reading columns from right to left, bottom to top
n = len(matrix)
rotated = []
for i in range(n):  # for each column from left to right
    new_row = []
    for j in range(n-1, -1, -1):  # for each row from bottom to top
        new_row.append(matrix[j][i])
    rotated.append(new_row)

# Convert to output format
result = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{result}>>>")