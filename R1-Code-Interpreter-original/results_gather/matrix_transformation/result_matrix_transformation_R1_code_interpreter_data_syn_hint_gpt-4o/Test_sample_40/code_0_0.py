# Original matrix
matrix = [
    ["S", "Z", "M"],
    ["S", "D", "F"],
    ["V", "W", "X"]
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = ','.join([','.join(row) for row in rotated_matrix])
print(result)