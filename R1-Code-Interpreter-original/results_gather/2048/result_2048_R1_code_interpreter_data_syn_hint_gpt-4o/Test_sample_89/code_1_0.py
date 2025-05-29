# Initial matrix
matrix = [
    [0, 2, 0, 4, 0],
    [0, 16, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [32, 0, 8, 0, 0],
    [0, 0, 16, 0, 0]
]

# Starting position
position = (1, 3)

# Movement sequence
moves = 'rurlrr'

# Define movement directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Simulate the movement
for move in moves:
    dx, dy = directions[move]
    new_x = position[0] + dx
    new_y = position[1] + dy
    # Check bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        # Check if the current number can combine with the number at the new position
        if matrix[new_x][new_y] == matrix[position[0]][position[1]]:
            # Combine the numbers
            matrix[new_x][new_y] *= 2
            # Set the current position to 0
            matrix[position[0]][position[1]] = 0
        elif matrix[position[0]][position[1]] != 0:
            # Move the number to the new position if it's not zero
            matrix[new_x][new_y] = matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
        # Update the position
        position = (new_x, new_y)

# Print the final matrix
print(matrix)