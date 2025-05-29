# Initial matrix
matrix = [
    [0, 16],
    [4, 2]
]

# Move sequence
moves = 'dl'

# Function to perform the moves
def move_number(matrix, moves):
    # Initial position of the number to move
    pos = (1, 1)  # (row, column) in 0-indexed
    for move in moves:
        if move == 'd':
            # Move down
            if pos[0] < 1:  # Check if not at the bottom row
                matrix[pos[0]][pos[1]], matrix[pos[0] + 1][pos[1]] = matrix[pos[0] + 1][pos[1]], matrix[pos[0]][pos[1]]
                pos = (pos[0] + 1, pos[1])
        elif move == 'l':
            # Move left
            if pos[1] > 0:  # Check if not at the leftmost column
                matrix[pos[0]][pos[1]], matrix[pos[0]][pos[1] - 1] = matrix[pos[0]][pos[1] - 1], matrix[pos[0]][pos[1]]
                pos = (pos[0], pos[1] - 1)
    return matrix

# Perform the moves
result_matrix = move_number(matrix, moves)

# Print the result
print(result_matrix)