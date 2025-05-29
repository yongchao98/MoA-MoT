# Input matrix as a string
matrix_str = """8 1 9 9 3 6 7 4 5
7 8 5 2 8 1 0 8 7
1 5 0 1 9 2 3 3 0
8 6 7 8 6 8 2 8 0
9 2 9 9 5 2 1 0 7
6 9 4 1 3 6 6 6 6
1 4 6 2 0 0 6 6 2
3 1 5 5 4 9 9 3 7
6 7 7 8 8 7 6 4 0"""

# Split into rows and process each row
matrix = [row.split() for row in matrix_str.strip().split('\n')]

# Format the output
output = '\n'.join(' '.join(row) for row in matrix)

print(output)