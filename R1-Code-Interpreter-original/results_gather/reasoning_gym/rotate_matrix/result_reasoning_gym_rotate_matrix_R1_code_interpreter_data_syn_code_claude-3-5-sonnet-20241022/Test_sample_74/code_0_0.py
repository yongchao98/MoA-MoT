# Parse the input matrix
input_str = """0 8 8 7 0 4 1 3 7
5 6 7 4 3 7 4 7 2
2 4 5 7 5 1 4 5 5
2 2 6 1 0 1 5 1 2
5 7 8 4 9 4 7 7 6
3 2 2 7 3 2 3 3 5
9 8 3 2 5 4 2 9 5
8 7 8 6 3 6 3 8 5
6 6 2 2 9 0 5 8 8"""

# Convert input string to matrix
matrix = [list(map(int, row.split())) for row in input_str.strip().split('\n')]

# Create rotated matrix
n = len(matrix)
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Perform rotation
for i in range(n):
    for j in range(n):
        rotated[j][n-1-i] = matrix[i][j]

# Format output
result = '\n'.join(' '.join(str(num) for num in row) for row in rotated)
print(f"<<<{result}>>>")