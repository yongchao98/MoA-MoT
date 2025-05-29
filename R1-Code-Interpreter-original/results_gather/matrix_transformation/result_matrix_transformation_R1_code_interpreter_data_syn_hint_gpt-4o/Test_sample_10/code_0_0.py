# Original matrix
matrix = [
    ['M', 'X'],
    ['B', 'U']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = [
    [matrix[1][0], matrix[0][0]],
    [matrix[1][1], matrix[0][1]]
]

# Flatten the rotated matrix to a list of values and symbols
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(flattened_result)