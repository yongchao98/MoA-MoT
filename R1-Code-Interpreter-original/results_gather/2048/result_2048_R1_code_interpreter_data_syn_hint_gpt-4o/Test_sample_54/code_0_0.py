# Initial matrix
matrix = [
    [4, 2, 0, 0],
    [0, 0, 0, 0],
    [2, 8, 0, 0],
    [32, 32, 4, 0]
]

# Initial position of the number
position = (3, 3)  # 0-indexed (4, 4) in 1-indexed

# Move sequence
moves = 'dluu'

# Perform the moves
for move in moves:
    if move == 'd':
        if position[0] < 3:  # Can move down
            position = (position[0] + 1, position[1])
    elif move == 'l':
        if position[1] > 0:  # Can move left
            position = (position[0], position[1] - 1)
    elif move == 'u':
        if position[0] > 0:  # Can move up
            position = (position[0] - 1, position[1])
    elif move == 'r':
        if position[1] < 3:  # Can move right
            position = (position[0], position[1] + 1)

# Since the number was 0, the matrix remains unchanged
print(matrix)