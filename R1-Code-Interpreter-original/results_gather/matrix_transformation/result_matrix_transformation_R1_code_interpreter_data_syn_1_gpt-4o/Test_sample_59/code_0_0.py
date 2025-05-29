# Original matrix
matrix = [
    ['Z', 'P', 'E'],
    ['Z', 'T', 'V'],
    ['G', 'I', 'V']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = []
n = len(matrix)
for i in range(n):
    new_row = [matrix[j][i] for j in range(n-1, -1, -1)]
    rotated_matrix.append(new_row)

# Flatten the rotated matrix to a single list
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)