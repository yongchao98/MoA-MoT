# Original matrix
matrix = [
    ['K', 'P', 'E', '3'],
    ['Y', '2', 'P', 'S'],
    ['L', 'J', 'U', '3']
]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse each row of the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Flatten the rotated matrix to match the required output format
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(flattened_result)