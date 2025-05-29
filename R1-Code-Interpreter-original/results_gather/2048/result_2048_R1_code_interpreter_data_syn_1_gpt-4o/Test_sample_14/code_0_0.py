# Initial matrix
matrix = [
    [0, 8],
    [2, 2]
]

# Move sequence
sequence = 'dl'

# Function to perform the move
def move_and_combine(matrix, sequence):
    # Start position of the number to move
    pos = (1, 1)  # (row, column) in zero-indexed

    for move in sequence:
        if move == 'd':
            # Move down
            if pos[0] < 1:  # Check if not at the bottom row
                pos = (pos[0] + 1, pos[1])
        elif move == 'l':
            # Move left
            if pos[1] > 0:  # Check if not at the leftmost column
                pos = (pos[0], pos[1] - 1)

    # After moving, check for combination
    if matrix[pos[0]][pos[1]] == 2:
        matrix[pos[0]][pos[1]] *= 2  # Combine to form 4
        matrix[1][1] = 0  # Clear the original position

    return matrix

# Perform the move and combination
result_matrix = move_and_combine(matrix, sequence)

# Print the resulting matrix
print(result_matrix)