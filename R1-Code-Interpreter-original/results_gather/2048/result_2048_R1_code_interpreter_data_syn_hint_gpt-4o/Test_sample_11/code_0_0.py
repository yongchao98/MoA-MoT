# Initial matrix
matrix = [
    [4, 4],
    [8, 8]
]

# Move left ('l') - no change as it's already at the leftmost position
# Move right ('r')
if matrix[0][0] == matrix[0][1]:
    matrix[0][1] *= 2
    matrix[0][0] = 0

# Print the final matrix
print(matrix)