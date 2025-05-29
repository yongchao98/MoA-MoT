# Original matrix
matrix = [
    [9, 2, 2],
    [6, 7, 7],
    [7, 7, 4]
]

# Step 1: Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Step 2: Reverse the order of the rows in the transposed matrix
rotated_matrix = [list(row) for row in transposed_matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)