# Initial matrix
matrix = [
    [0, 16],
    [2, 0]
]

# Move left ('l')
matrix[1][0], matrix[1][1] = matrix[1][1], matrix[1][0]

# Move up ('u')
matrix[0][0], matrix[1][0] = matrix[1][0], matrix[0][0]

# Print the final matrix
print(matrix)