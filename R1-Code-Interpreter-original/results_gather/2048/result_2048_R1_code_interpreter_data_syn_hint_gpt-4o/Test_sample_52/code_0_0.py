# Initial matrix
matrix = [
    [0, 2, 0, 0],
    [0, 0, 32, 2],
    [0, 0, 4, 0],
    [0, 0, 0, 0]
]

# Initial position of the number
position = (1, 0)  # (row, column) in 0-indexed

# Movement sequence
sequence = 'dddl'

# Function to move the number
def move_number(matrix, position, sequence):
    for move in sequence:
        row, col = position
        if move == 'd' and row < 3:
            position = (row + 1, col)
        elif move == 'u' and row > 0:
            position = (row - 1, col)
        elif move == 'l' and col > 0:
            position = (row, col - 1)
        elif move == 'r' and col < 3:
            position = (row, col + 1)
    return matrix

# Apply the movement
final_matrix = move_number(matrix, position, sequence)

# Print the final matrix
print(final_matrix)