# Define the original matrix
matrix = [
    ['7', '1', '1'],
    ['5', '3', '8']
]

# Rotate the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Convert tuples to lists and join elements to form the final output
result = [','.join(row) for row in rotated_matrix]

# Print the result
print(result)