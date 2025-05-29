# Original matrix
matrix = [
    ['8', 'W', 'N', 'K'],
    ['P', 'U', 'F', '7'],
    ['R', 'L', '3', 'F']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = [','.join(row) for row in rotated_matrix]
print(result)