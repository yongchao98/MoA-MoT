# Initial matrix
matrix = [
    [0, 32, 0, 4, 0],
    [0, 0, 0, 4, 0],
    [0, 0, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [16, 32, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 4)  # (5, 5) in 0-indexed

# Movement sequence
moves = 'uldll'

# Movement directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Simulate the moves
for move in moves:
    dx, dy = directions[move]
    new_position = (position[0] + dx, position[1] + dy)
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        position = new_position

# The number being moved is 0, so no change occurs in the matrix
# Print the final matrix
for row in matrix:
    print(row)