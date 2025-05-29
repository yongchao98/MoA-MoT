# Original matrix
matrix = [
    ['C', 'Z'],
    ['I', 'U']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output as a list of values and symbols
result = [item for sublist in rotated_matrix for item in sublist]
print(result)