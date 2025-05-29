# Initial matrix
matrix = [
    [0, 0, 4],
    [16, 4, 0],
    [0, 0, 0]
]

# Initial position of the number to move
position = (2, 1)  # (3, 2) in 0-indexed

# Sequence of moves
moves = 'rll'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'r':
            if y < 2:  # Move right
                y += 1
        elif move == 'l':
            if y > 0:  # Move left
                y -= 1
        elif move == 'u':
            if x > 0:  # Move up
                x -= 1
        elif move == 'd':
            if x < 2:  # Move down
                x += 1
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)