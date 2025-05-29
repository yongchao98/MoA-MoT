# Initial matrix
matrix = [
    [8, 0, 0, 32, 0],
    [4, 0, 0, 16, 0],
    [0, 0, 0, 0, 0],
    [2, 4, 0, 2, 4],
    [4, 0, 0, 0, 0]
]

# Initial position of the number to move
position = (2, 4)

# Sequence of moves
moves = 'dldur'

# Function to move the position
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
    position = move(position, move_direction)

# Print the final matrix
print(matrix)