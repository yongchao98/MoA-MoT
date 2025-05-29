# Initial matrix
matrix = [
    [2, 0, 0],
    [0, 2, 0],
    [0, 16, 16]
]

# Initial position of the number to move
position = (2, 0)  # (row, column) in 0-indexed

# Move sequence
moves = 'dld'

# Function to move the number
def move_number(matrix, position, moves):
    for move in moves:
        if move == 'd':
            if position[0] < 2:  # Can move down if not in the last row
                position = (position[0] + 1, position[1])
        elif move == 'l':
            if position[1] > 0:  # Can move left if not in the first column
                position = (position[0], position[1] - 1)
        elif move == 'r':
            if position[1] < 2:  # Can move right if not in the last column
                position = (position[0], position[1] + 1)
        elif move == 'u':
            if position[0] > 0:  # Can move up if not in the first row
                position = (position[0] - 1, position[1])
    return matrix

# Apply the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)