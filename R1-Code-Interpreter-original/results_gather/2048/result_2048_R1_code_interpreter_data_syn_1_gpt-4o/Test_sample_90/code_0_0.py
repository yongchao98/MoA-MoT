# Initial grid
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 16],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [4, 2, 0, 16, 32]
]

# Initial position of the number to move
position = (2, 3)

# Move sequence
moves = 'ldullu'

# Move deltas
move_deltas = {
    'l': (0, -1),
    'r': (0, 1),
    'u': (-1, 0),
    'd': (1, 0)
}

# Simulate the moves
for move in moves:
    delta = move_deltas[move]
    new_position = (position[0] + delta[0], position[1] + delta[1])
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        position = new_position

# Since the initial number is 0, the grid remains unchanged
print(grid)