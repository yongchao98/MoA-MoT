# Initial matrix
matrix = [
    [0, 0, 0, 0],
    [0, 0, 16, 32],
    [16, 32, 4, 4],
    [0, 4, 0, 0]
]

# Initial position of the number to move
position = (3, 1)

# Sequence of moves
moves = 'uuud'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'u':
            x = max(0, x - 1)
        elif move == 'd':
            x = min(3, x + 1)
        elif move == 'l':
            y = max(0, y - 1)
        elif move == 'r':
            y = min(3, y + 1)
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)