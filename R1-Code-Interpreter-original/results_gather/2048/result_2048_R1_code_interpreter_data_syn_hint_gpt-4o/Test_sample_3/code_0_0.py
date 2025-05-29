# Initial matrix
matrix = [
    [32, 2],
    [32, 4]
]

# Move 'u' (up)
matrix[0][1], matrix[1][1] = matrix[1][1], matrix[0][1]

# Move 'r' (right) - no change since it's already at the rightmost column

# Print the final matrix
print(matrix)