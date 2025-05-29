# Initial matrix
matrix = [
    [0, 4, 0],
    [8, 0, 2],
    [8, 2, 0]
]

# Initial position of the number to move
position = (1, 2)

# Move sequence
moves = 'rdl'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'r':  # Move right
            if y < 2:
                y += 1
        elif move == 'd':  # Move down
            if x < 2:
                x += 1
        elif move == 'l':  # Move left
            if y > 0:
                y -= 1
        elif move == 'u':  # Move up
            if x > 0:
                x -= 1

        # Check for combination
        if matrix[x][y] == 0:
            matrix[x][y] = matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
        elif matrix[x][y] == matrix[position[0]][position[1]]:
            matrix[x][y] *= 2
            matrix[position[0]][position[1]] = 0
        else:
            # If the move is not possible, revert the position
            if move == 'r':
                y -= 1
            elif move == 'd':
                x -= 1
            elif move == 'l':
                y += 1
            elif move == 'u':
                x += 1

        # Update position
        position = (x, y)

    return matrix

# Perform the move
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)