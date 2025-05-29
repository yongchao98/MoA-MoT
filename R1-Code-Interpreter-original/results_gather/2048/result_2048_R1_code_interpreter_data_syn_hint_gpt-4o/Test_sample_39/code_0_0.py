# Initial matrix
matrix = [
    [8, 4, 4],
    [0, 16, 0],
    [16, 0, 0]
]

# Initial position of the number to move
position = (1, 1)

# Move sequence
moves = 'llr'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'l' and y > 0:
            y -= 1
        elif move == 'r' and y < 2:
            y += 1
        elif move == 'u' and x > 0:
            x -= 1
        elif move == 'd' and x < 2:
            x += 1
    return x, y

# Perform the moves
new_position = move_number(matrix, position, moves)

# Update the matrix
matrix[new_position[0]][new_position[1]] = matrix[position[0]][position[1]]
matrix[position[0]][position[1]] = 0

# Print the updated matrix
print(matrix)