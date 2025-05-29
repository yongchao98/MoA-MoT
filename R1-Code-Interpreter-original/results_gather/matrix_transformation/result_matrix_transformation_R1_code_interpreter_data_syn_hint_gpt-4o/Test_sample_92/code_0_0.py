# Original matrix
matrix = [
    ['*', '#', '$', '*'],
    ['%', '@', '#', '@'],
    ['&', '*', '^', '%'],
    ['^', '!', '$', '*']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Flatten the rotated matrix for the output format
flattened_rotated_matrix = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(flattened_rotated_matrix)