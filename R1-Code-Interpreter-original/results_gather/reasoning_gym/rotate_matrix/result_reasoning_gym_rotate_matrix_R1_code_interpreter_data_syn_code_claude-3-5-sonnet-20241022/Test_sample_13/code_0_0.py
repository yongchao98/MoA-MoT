# Input matrix as a string
matrix_str = """6 7 4 1 4 5
0 1 7 0 1 4
1 4 0 9 1 2
4 2 3 6 1 9
6 2 0 0 2 5
4 5 5 1 6 1"""

# Split into rows and convert to integers
matrix = [list(map(int, row.split())) for row in matrix_str.splitlines()]

# For 0-degree rotation, matrix stays the same
# Format the output
output = "\n".join(" ".join(str(num) for num in row) for row in matrix)

print(f"<<<{output}>>>")