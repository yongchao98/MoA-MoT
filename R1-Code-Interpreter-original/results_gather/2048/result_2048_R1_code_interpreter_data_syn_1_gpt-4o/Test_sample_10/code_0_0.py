# Initial matrix
matrix = [
    [8, 4],
    [16, 0]
]

# Move sequence
moves = 'lr'

# Function to perform the moves
def move_number(matrix, moves):
    # Start position of the number to move
    pos = (1, 1)  # (row, column) in zero-indexed
    for move in moves:
        if move == 'l':  # Move left
            if pos[1] > 0:
                # Swap the number with the left cell
                matrix[pos[0]][pos[1]], matrix[pos[0]][pos[1] - 1] = matrix[pos[0]][pos[1] - 1], matrix[pos[0]][pos[1]]
                pos = (pos[0], pos[1] - 1)
        elif move == 'r':  # Move right
            if pos[1] < 1:
                # Swap the number with the right cell
                matrix[pos[0]][pos[1]], matrix[pos[0]][pos[1] + 1] = matrix[pos[0]][pos[1] + 1], matrix[pos[0]][pos[1]]
                pos = (pos[0], pos[1] + 1)
    return matrix

# Perform the moves
result_matrix = move_number(matrix, moves)

# Print the resulting matrix
print(result_matrix)