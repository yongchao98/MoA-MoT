# Initial matrix
matrix = [
    [0, 0, 4, 0],
    [0, 0, 0, 0],
    [0, 0, 4, 0],
    [2, 0, 16, 16]
]

# Move sequence: 'urdu'
# Initial position of the number is (4, 4) with value 16
position = (3, 3)  # zero-indexed position (4, 4)

# Move 'u' (up)
if position[0] > 0:
    if matrix[position[0] - 1][position[1]] == matrix[position[0]][position[1]]:
        matrix[position[0] - 1][position[1]] *= 2
        matrix[position[0]][position[1]] = 0
    else:
        matrix[position[0] - 1][position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
    position = (position[0] - 1, position[1])

# Move 'r' (right)
# Since the number is already at the rightmost column, this move does not change the position

# Move 'd' (down)
if position[0] < 3:
    if matrix[position[0] + 1][position[1]] == matrix[position[0]][position[1]]:
        matrix[position[0] + 1][position[1]] *= 2
        matrix[position[0]][position[1]] = 0
    else:
        matrix[position[0] + 1][position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
    position = (position[0] + 1, position[1])

# Move 'u' (up)
if position[0] > 0:
    if matrix[position[0] - 1][position[1]] == matrix[position[0]][position[1]]:
        matrix[position[0] - 1][position[1]] *= 2
        matrix[position[0]][position[1]] = 0
    else:
        matrix[position[0] - 1][position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
    position = (position[0] - 1, position[1])

# Print the final matrix
print(matrix)