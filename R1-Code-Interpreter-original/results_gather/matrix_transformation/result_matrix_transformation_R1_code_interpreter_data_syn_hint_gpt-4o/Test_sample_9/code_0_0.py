# Original matrix
matrix = [
    ['S', 'V'],
    ['U', 'O']
]

# Step 1: Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Step 2: Reverse each row in the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Print the result in the required format
print(rotated_matrix)