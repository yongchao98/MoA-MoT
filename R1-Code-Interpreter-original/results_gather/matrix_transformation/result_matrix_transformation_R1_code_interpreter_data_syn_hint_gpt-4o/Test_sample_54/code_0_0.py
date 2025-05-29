# Original matrix
matrix = [
    ['W', 'V', 'Z'],
    ['P', 'P', 'X'],
    ['O', 'J', 'V']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = ','.join([','.join(row) for row in rotated_matrix])
print(result)