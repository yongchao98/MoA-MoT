# Initial matrix
matrix = [
    [32, 2],
    [32, 4]
]

# Move sequence
moves = 'ur'

# Function to perform the moves
def move_2048(matrix, moves):
    # Start with the position of the number 4
    pos = (1, 1)  # (row, column) in zero-indexed

    for move in moves:
        if move == 'u':
            # Move up
            if pos[0] > 0:
                # Swap the numbers
                matrix[pos[0]][pos[1]], matrix[pos[0] - 1][pos[1]] = matrix[pos[0] - 1][pos[1]], matrix[pos[0]][pos[1]]
                pos = (pos[0] - 1, pos[1])
        elif move == 'r':
            # Move right
            if pos[1] < len(matrix[0]) - 1:
                # Swap the numbers
                matrix[pos[0]][pos[1]], matrix[pos[0]][pos[1] + 1] = matrix[pos[0]][pos[1] + 1], matrix[pos[0]][pos[1]]
                pos = (pos[0], pos[1] + 1)

    return matrix

# Perform the moves
result_matrix = move_2048(matrix, moves)

# Print the result
print(result_matrix)