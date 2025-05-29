# Original matrix
matrix = [
    ['Z', 'C'],
    ['R', 'U']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = [item for sublist in rotated_matrix for item in sublist]
print(result)