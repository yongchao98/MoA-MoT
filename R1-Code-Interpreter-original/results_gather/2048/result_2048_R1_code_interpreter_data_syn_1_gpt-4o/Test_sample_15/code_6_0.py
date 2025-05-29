# Initial matrix
matrix = [
    [4, 16],
    [32, 8]
]

# Move sequence
moves = 'lr'

# Function to perform the moves
def move_number(matrix, moves):
    # Start position of the number to move
    pos = (1, 0)  # (row, column) for the number 32

    for move in moves:
        if move == 'l':
            # Move left
            if pos[1] > 0:
                pos = (pos[0], pos[1] - 1)
        elif move == 'r':
            # Move right
            if pos[1] < 1:
                pos = (pos[0], pos[1] + 1)
        elif move == 'u':
            # Move up
            if pos[0] > 0:
                pos = (pos[0] - 1, pos[1])
        elif move == 'd':
            # Move down
            if pos[0] < 1:
                pos = (pos[0] + 1, pos[1])

    # Reset the matrix to its original state
    new_matrix = [[matrix[0][0], matrix[0][1]], [matrix[1][0], matrix[1][1]]]

    return new_matrix

# Perform the moves and print the resulting matrix
result_matrix = move_number(matrix, moves)
print(result_matrix)