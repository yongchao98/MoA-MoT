# Original matrix
matrix = [
    ["M", "J", "B"],
    ["T", "V", "R"],
    ["M", "B", "E"]
]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse each row of the transposed matrix to get the rotated matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Formatting the output
result = ','.join([','.join(row) for row in rotated_matrix])
print(result)