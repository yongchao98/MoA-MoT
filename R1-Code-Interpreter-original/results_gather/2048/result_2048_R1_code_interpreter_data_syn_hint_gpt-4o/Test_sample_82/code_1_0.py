# Initial matrix
matrix = [
    [0, 8, 16, 2, 2],
    [2, 8, 0, 0, 0],
    [0, 4, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [8, 0, 8, 0, 0]
]

# Initial position of the number to move
position = (1, 3)

# Sequence of moves
moves = 'uldlur'

# Function to perform a move
def move(position, direction):
    x, y = position
    if direction == 'u':
        return (x - 1, y)
    elif direction == 'd':
        return (x + 1, y)
    elif direction == 'l':
        return (x, y - 1)
    elif direction == 'r':
        return (x, y + 1)
    return position

# Simulate the moves
for move_direction in moves:
    new_position = move(position, move_direction)
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        # Move the number
        number = matrix[position[0]][position[1]]
        target_number = matrix[new_position[0]][new_position[1]]
        # Check if they can combine
        if number != 0 and number == target_number:
            matrix[new_position[0]][new_position[1]] *= 2
            matrix[position[0]][position[1]] = 0
        elif number != 0:
            # Just move the number
            matrix[new_position[0]][new_position[1]] = number
            matrix[position[0]][position[1]] = 0
        # Update the position
        position = new_position

# Print the final matrix
print(matrix)