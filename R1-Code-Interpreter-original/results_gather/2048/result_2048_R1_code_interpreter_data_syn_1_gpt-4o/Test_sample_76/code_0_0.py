# Initial grid
grid = [
    [0, 32, 0, 4, 0],
    [0, 0, 0, 4, 0],
    [0, 0, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [16, 32, 0, 0, 0]
]

# Initial position of the number
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
    position = new_position

# Print the final grid
for row in grid:
    print(row)