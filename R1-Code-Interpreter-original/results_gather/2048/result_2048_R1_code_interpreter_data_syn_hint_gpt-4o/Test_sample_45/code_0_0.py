# Initial matrix
matrix = [
    [0, 0, 0, 0],
    [0, 0, 16, 32],
    [16, 32, 4, 4],
    [0, 4, 0, 0]
]

# Sequence of moves
moves = 'uuud'

# Initial position of the number to move
position = [3, 1]

# Execute the moves
for move in moves:
    if move == 'u':
        position[0] -= 1
    elif move == 'd':
        position[0] += 1
    elif move == 'l':
        position[1] -= 1
    elif move == 'r':
        position[1] += 1

# The number is 0, so it doesn't change the matrix
# Print the final matrix
print(matrix)