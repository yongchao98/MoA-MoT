# Original matrix
matrix = [
    ['L', 'G', 'J', 'M'],
    ['B', 'I', 'F', 'K'],
    ['F', '7', 'R', '3']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Flattening the rotated matrix and formatting the output
result = ','.join([element for row in rotated_matrix for element in row])
print(result)