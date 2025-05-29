# Original matrix
matrix = [
    ['3', 'X', 'S', 'N'],
    ['3', 'D', 'L', 'O'],
    ['B', '6', 'C', '3']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Flatten the rotated matrix to a single list
flattened_rotated_matrix = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(flattened_rotated_matrix)