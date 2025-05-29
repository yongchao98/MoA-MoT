# Initial matrix
matrix = [
    [0, 32, 8, 0],
    [0, 0, 0, 0],
    [8, 0, 0, 32],
    [16, 0, 0, 0]
]

# Initial position of the number
position = (2, 1)

# Movement sequence
moves = 'dlul'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'd':  # Move down
            if x < 3 and matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
                x += 1
        elif move == 'l':  # Move left
            if y > 0 and matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
                y -= 1
        elif move == 'u':  # Move up
            if x > 0 and matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
                x -= 1
        elif move == 'r':  # Move right
            if y < 3 and matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
                y += 1
    return matrix

# Perform the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)