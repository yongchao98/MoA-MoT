# Original matrix
matrix = [
    ['L', 'G', 'J', 'M'],
    ['B', 'I', 'F', 'K'],
    ['F', '7', 'R', '3']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = ','.join([','.join(row) for row in rotated_matrix])
print(result)