# Initial matrix
matrix = [
    [0, 0, 0],
    [16, 16, 0],
    [16, 0, 0]
]

# Initial position of the number to move
position = (2, 2)

# Move sequence
moves = 'rud'

# Function to perform a move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'r':  # Move right
        if y < 2:
            y += 1
    elif direction == 'u':  # Move up
        if x > 0:
            x -= 1
    elif direction == 'd':  # Move down
        if x < 2:
            x += 1
    elif direction == 'l':  # Move left
        if y > 0:
            y -= 1
    return (x, y)

# Simulate the moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)