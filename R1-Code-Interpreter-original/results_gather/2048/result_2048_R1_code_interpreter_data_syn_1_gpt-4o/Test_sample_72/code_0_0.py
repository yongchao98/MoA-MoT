# Initial matrix
matrix = [
    [0, 0, 4, 0, 0],
    [0, 0, 32, 0, 32],
    [0, 0, 0, 0, 0],
    [8, 0, 0, 0, 16],
    [0, 32, 0, 0, 4]
]

# Initial position of the number to move
position = (4, 3)

# Move sequence
moves = 'lluud'

# Function to perform a move
def move(position, direction):
    x, y = position
    if direction == 'l':
        return (x, max(0, y - 1))
    elif direction == 'r':
        return (x, min(4, y + 1))
    elif direction == 'u':
        return (max(0, x - 1), y)
    elif direction == 'd':
        return (min(4, x + 1), y)

# Simulate the moves
for move_direction in moves:
    position = move(position, move_direction)

# Print the final matrix
print(matrix)