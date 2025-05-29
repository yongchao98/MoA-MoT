# Initial matrix
matrix = [
    [0, 2, 0, 0],
    [0, 0, 32, 2],
    [0, 0, 4, 0],
    [0, 0, 0, 0]
]

# Initial position of the number to move
position = (1, 0)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'dddl'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'd' and row < 3:
            position = (row + 1, col)
        elif move == 'u' and row > 0:
            position = (row - 1, col)
        elif move == 'l' and col > 0:
            position = (row, col - 1)
        elif move == 'r' and col < 3:
            position = (row, col + 1)
    return position

# Get the final position after all moves
final_position = move_number(matrix, position, moves)

# Move the number in the matrix
matrix[final_position[0]][final_position[1]] = matrix[position[0]][position[1]]
matrix[position[0]][position[1]] = 0

# Print the final matrix
print(matrix)