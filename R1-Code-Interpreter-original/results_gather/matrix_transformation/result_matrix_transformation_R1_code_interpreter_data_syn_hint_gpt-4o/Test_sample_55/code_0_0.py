# Original matrix
matrix = [
    ['R', 'G', 'J'],
    ['Z', 'U', 'L'],
    ['K', 'S', 'V']
]

# Rotate the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Convert tuples to lists and format the output
rotated_matrix = [list(row) for row in rotated_matrix]
formatted_output = ','.join([','.join(row) for row in rotated_matrix])

print(formatted_output)