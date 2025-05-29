# Initial matrix
matrix = [
    [32, 0, 0, 0, 0],
    [0, 4, 4, 2, 0],
    [0, 2, 0, 0, 0],
    [0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 2)
value = matrix[position[0]][position[1]]

# Sequence of moves
moves = 'lrdlr'

# Function to update position based on move
def update_position(pos, move):
    if move == 'l':
        return (pos[0], max(0, pos[1] - 1))
    elif move == 'r':
        return (pos[0], min(4, pos[1] + 1))
    elif move == 'u':
        return (max(0, pos[0] - 1), pos[1])
    elif move == 'd':
        return (min(4, pos[0] + 1), pos[1])

# Simulate the moves
for move in moves:
    new_position = update_position(position, move)
    # Move the value to the new position
    matrix[position[0]][position[1]] = 0
    matrix[new_position[0]][new_position[1]] = value
    position = new_position

# Print the final matrix
print(matrix)