# Original matrix
matrix = [
    ['F', 'J'],
    ['U', 'T']
]

# Step 1: Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Step 2: Reverse each row in the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Prepare the result in the required format
result = ','.join([','.join(row) for row in rotated_matrix])
print(result)