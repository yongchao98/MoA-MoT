# Initial matrix
matrix = [
    [0, 4, 0, 4],
    [0, 0, 0, 0],
    [0, 4, 16, 0],
    [0, 0, 0, 0]
]

# Position of the number to move
position = (2, 3)  # (row, column) in 0-indexed

# Move sequence
moves = 'duuu'

# Function to move the number
def move_number(matrix, position, moves):
    row, col = position
    for move in moves:
        if move == 'd' and row < 3:
            row += 1
        elif move == 'u' and row > 0:
            row -= 1
        elif move == 'l' and col > 0:
            col -= 1
        elif move == 'r' and col < 3:
            col += 1
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the result
print(result_matrix)