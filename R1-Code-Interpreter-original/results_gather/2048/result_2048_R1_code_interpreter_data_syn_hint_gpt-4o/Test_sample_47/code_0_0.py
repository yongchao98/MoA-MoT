# Initial matrix
matrix = [
    [16, 8, 0, 2],
    [0, 0, 32, 0],
    [0, 4, 32, 0],
    [0, 32, 0, 2]
]

# Initial position of the number to move
position = (3, 3)  # (row, column) in 0-indexed

# Move sequence
moves = 'ddud'

# Function to move the number according to the sequence
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'u' and row > 0:
            position = (row - 1, col)
        elif move == 'd' and row < 3:
            position = (row + 1, col)
        elif move == 'l' and col > 0:
            position = (row, col - 1)
        elif move == 'r' and col < 3:
            position = (row, col + 1)
    return position

# Get the final position
final_position = move_number(matrix, position, moves)

# Update the matrix
final_matrix = [row[:] for row in matrix]  # Copy the matrix
final_matrix[3][3] = 0  # Clear the initial position
final_matrix[final_position[0]][final_position[1]] = 2  # Place the number in the final position

print(final_matrix)