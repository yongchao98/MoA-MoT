# Original matrix
matrix = [
    [8, 8, 1, 2, 6, 3, 4, 0],
    [3, 1, 9, 0, 1, 2, 8, 4],
    [6, 9, 6, 5, 5, 1, 5, 4],
    [9, 2, 1, 8, 1, 9, 1, 4],
    [5, 1, 4, 0, 5, 6, 1, 7],
    [7, 3, 3, 2, 4, 3, 0, 0],
    [6, 0, 5, 5, 7, 7, 9, 8],
    [2, 3, 7, 7, 5, 9, 0, 4]
]

# Step 1: Reverse the order of the rows
matrix.reverse()

# Step 2: Reverse the order of elements in each row
for row in matrix:
    row.reverse()

# Print the rotated matrix
print(matrix)