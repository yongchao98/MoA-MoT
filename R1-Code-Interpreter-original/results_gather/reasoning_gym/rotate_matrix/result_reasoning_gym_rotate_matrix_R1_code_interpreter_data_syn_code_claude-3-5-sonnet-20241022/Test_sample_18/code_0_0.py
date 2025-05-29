# Input matrix as a string
matrix_str = """4 5 3 7 5 1 8 5 3
8 1 2 5 1 5 5 3 3
1 5 9 9 6 1 1 9 5
6 0 3 5 4 6 1 3 7
4 9 7 2 4 9 3 8 9
1 6 9 1 8 7 6 8 8
2 7 9 8 4 1 5 8 2
0 5 0 9 6 7 8 3 2
9 1 2 4 6 3 5 0 8"""

# Split into rows and convert to matrix
matrix = [row.split() for row in matrix_str.strip().split('\n')]

# Format output
output = '\n'.join(' '.join(row) for row in matrix)
print(output)