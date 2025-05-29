# Initial matrix
matrix = [
    [2, 8, 0],
    [4, 0, 0],
    [0, 8, 16]
]

# Initial position of the number to move
position = (1, 2)

# Sequence of moves
moves = 'url'

# Function to perform the move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'u':
        if x > 0:
            x -= 1
    elif direction == 'd':
        if x < 2:
            x += 1
    elif direction == 'l':
        if y > 0:
            y -= 1
    elif direction == 'r':
        if y < 2:
            y += 1
    return (x, y)

# Simulate the moves
for move in moves:
    new_position = move_number(matrix, position, move)
    # Move the number
    matrix[new_position[0]][new_position[1]] = matrix[position[0]][position[1]]
    matrix[position[0]][position[1]] = 0
    position = new_position

print(matrix)