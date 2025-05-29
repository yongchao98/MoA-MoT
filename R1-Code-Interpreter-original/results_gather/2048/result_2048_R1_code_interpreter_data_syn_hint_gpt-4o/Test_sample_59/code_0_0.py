# Initial matrix
matrix = [
    [32, 8, 32, 0],
    [0, 0, 32, 0],
    [0, 8, 0, 0],
    [4, 0, 0, 0]
]

# Initial position of the number
position = (2, 2)

# Move sequence
moves = 'rudl'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'r':  # Move right
            if y + 1 < 4 and (matrix[x][y + 1] == 0 or matrix[x][y + 1] == matrix[x][y]):
                y += 1
        elif move == 'u':  # Move up
            if x - 1 >= 0 and (matrix[x - 1][y] == 0 or matrix[x - 1][y] == matrix[x][y]):
                x -= 1
        elif move == 'd':  # Move down
            if x + 1 < 4 and (matrix[x + 1][y] == 0 or matrix[x + 1][y] == matrix[x][y]):
                x += 1
        elif move == 'l':  # Move left
            if y - 1 >= 0 and (matrix[x][y - 1] == 0 or matrix[x][y - 1] == matrix[x][y]):
                y -= 1

    # Update the matrix with the final position of the number
    matrix[position[0]][position[1]] = 0  # Clear the initial position
    matrix[x][y] = 8  # Place the number in the final position

    return matrix

# Perform the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)