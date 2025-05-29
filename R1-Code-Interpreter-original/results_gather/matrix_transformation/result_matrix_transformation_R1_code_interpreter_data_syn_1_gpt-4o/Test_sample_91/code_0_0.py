# Original matrix
matrix = [
    ['#', '^', '#', '#'],
    ['%', '*', '%', '$'],
    ['&', '!', '&', '*'],
    ['%', '&', '%', '!']
]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse each row of the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Flatten the rotated matrix to a single list
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)