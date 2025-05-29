# Initial matrix
matrix = [
    [0, 4, 0, 2, 0],
    [0, 16, 8, 0, 32],
    [0, 0, 0, 16, 0],
    [4, 16, 0, 0, 0],
    [0, 8, 0, 0, 0]
]

# Initial position
position = (4, 4)

# Movement sequence
moves = 'ulurrl'

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

# Simulate the movement
for move_direction in moves:
    new_position = move(position, move_direction)
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        position = new_position

# Since the initial number was 0, the matrix remains unchanged
print(matrix)