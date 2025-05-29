# Original matrix
matrix = [
    ['2', 'J', 'C', '3'],
    ['E', '9', 'Y', 'A'],
    ['W', 'E', 'A', 'V']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = [','.join(row) for row in rotated_matrix]
print(result)