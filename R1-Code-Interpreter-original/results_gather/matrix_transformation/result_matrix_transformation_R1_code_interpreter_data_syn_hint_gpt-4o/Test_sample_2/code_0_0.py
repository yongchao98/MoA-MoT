# Original matrix
matrix = [
    ['W', 'U'],
    ['W', 'F']
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output as a list of strings
result = [','.join(row) for row in rotated_matrix]

# Printing the result
print(result)