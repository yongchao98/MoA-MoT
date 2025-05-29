# Original matrix
matrix = [
    ['^', '^', '@'],
    ['@', '%', '@', '&'],
    ['%', '*', '^', '%'],
    ['!', '!', '!', '*']
]

# Determine the maximum row length
max_length = max(len(row) for row in matrix)

# Pad shorter rows with None to make the matrix rectangular
padded_matrix = [row + [None] * (max_length - len(row)) for row in matrix]

# Transpose the padded matrix
transposed_matrix = list(zip(*padded_matrix))

# Reverse each row in the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Flatten the rotated matrix to a single list, ignoring None values
flattened_result = [item for sublist in rotated_matrix for item in sublist if item is not None]

# Print the result
print(flattened_result)