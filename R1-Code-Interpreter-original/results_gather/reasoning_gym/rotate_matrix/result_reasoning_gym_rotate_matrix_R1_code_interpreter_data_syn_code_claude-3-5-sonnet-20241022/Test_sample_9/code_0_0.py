# Input matrix
input_str = """5 3 2 6 9 4 8 4
6 9 0 3 0 7 2 2
5 2 4 9 9 0 1 9
0 6 4 3 1 5 7 2
8 0 1 1 2 2 0 4
5 3 9 0 9 9 2 6
0 4 6 2 1 1 6 0
7 6 3 0 1 2 1 7"""

# Convert input string to matrix
matrix = [list(map(int, row.split())) for row in input_str.splitlines()]

# For 180-degree rotation, we can simply reverse the matrix both horizontally and vertically
# or we can directly map each element to its opposite position
n = len(matrix)
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Rotate 180 degrees
for i in range(n):
    for j in range(n):
        rotated[i][j] = matrix[n-1-i][n-1-j]

# Format output
output = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{output}>>>")