# Initial matrix
matrix = [
    [2, 4, 32, 0],
    [0, 0, 16, 0],
    [0, 0, 0, 0],
    [2, 8, 0, 0]
]

# Move sequence
moves = 'lrud'

# Function to perform the moves
def move_number(matrix, moves):
    # Start position
    pos = (1, 3)  # (row, column) for (1, 4) in 0-indexed
    for move in moves:
        if move == 'l':
            pos = (pos[0], max(0, pos[1] - 1))
        elif move == 'r':
            pos = (pos[0], min(3, pos[1] + 1))
        elif move == 'u':
            pos = (max(0, pos[0] - 1), pos[1])
        elif move == 'd':
            pos = (min(3, pos[0] + 1), pos[1])
    return matrix

# Perform the moves
result_matrix = move_number(matrix, moves)

# Print the result
print(result_matrix)