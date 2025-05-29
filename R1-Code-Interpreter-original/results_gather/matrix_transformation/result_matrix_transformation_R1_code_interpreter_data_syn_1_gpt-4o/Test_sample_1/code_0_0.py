# Original matrix
matrix = [
    ['V', 'C'],
    ['Y', 'E']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Flatten the rotated matrix to a list of strings
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)