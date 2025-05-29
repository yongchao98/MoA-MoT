# Original matrix
matrix = [
    ['M', 'C', 'K', 'X'],
    ['8', '7', 'N', 'C'],
    ['Q', 'M', 'V', '1']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = ','.join([','.join(row) for row in rotated_matrix])
print(result)