# Initial matrix
matrix = [
    [0, 0, 0, 0, 8],
    [0, 0, 32, 0, 0],
    [0, 0, 0, 0, 0],
    [16, 0, 0, 0, 0],
    [0, 0, 0, 32, 8]
]

# Initial position of the number to move
position = (3, 1)

# Sequence of moves
moves = 'rlrrlr'

# Function to perform the move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'r':
        if y < 4:
            y += 1
    elif direction == 'l':
        if y > 0:
            y -= 1
    elif direction == 'u':
        if x > 0:
            x -= 1
    elif direction == 'd':
        if x < 4:
            x += 1
    return (x, y)

# Simulate the moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)