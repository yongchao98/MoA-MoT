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

# For 180-degree rotation:
# 1. Reverse the rows
# 2. Reverse each row
rotated = matrix[::-1]
for i in range(len(rotated)):
    rotated[i] = rotated[i][::-1]

# Format output
output = '\n'.join(' '.join(str(x) for x in row) for row in rotated)
print(f"<<<{output}>>>")