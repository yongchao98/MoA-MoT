# Original matrix as a list of strings
matrix = [
    "^ ^ @",
    "@ % @ &",
    "% * ^ %",
    "! ! ! *"
]

# Split each string into a list of characters
matrix = [row.split() for row in matrix]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse each row of the transposed matrix
rotated_matrix = [list(reversed(row)) for row in transposed_matrix]

# Flatten the rotated matrix into a single list
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)