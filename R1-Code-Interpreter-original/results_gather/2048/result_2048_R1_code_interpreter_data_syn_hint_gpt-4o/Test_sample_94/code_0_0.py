# Initial matrix
matrix = [
    [0, 0, 0, 0, 32],
    [0, 0, 0, 32, 32],
    [8, 0, 0, 0, 2],
    [2, 8, 4, 0, 0],
    [0, 4, 0, 0, 0]
]

# Initial position of the number to move
position = (2, 4)

# Sequence of moves
moves = 'ldldul'

# Function to perform a move
def move(position, direction):
    x, y = position
    if direction == 'l':
        y -= 1
    elif direction == 'r':
        y += 1
    elif direction == 'u':
        x -= 1
    elif direction == 'd':
        x += 1
    return x, y

# Perform the moves
for move_direction in moves:
    new_position = move(position, move_direction)
    x, y = position
    new_x, new_y = new_position
    # Move the number
    matrix[new_x][new_y] = matrix[x][y]
    matrix[x][y] = 0
    # Update the position
    position = new_position

# Print the final matrix
print(matrix)