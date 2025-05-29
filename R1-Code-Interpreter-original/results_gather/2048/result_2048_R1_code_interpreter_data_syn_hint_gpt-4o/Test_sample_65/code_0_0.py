# Initial matrix
matrix = [
    [0, 0, 0, 8, 0],
    [0, 16, 0, 16, 0],
    [0, 0, 0, 0, 8],
    [0, 2, 2, 0, 16],
    [2, 0, 0, 2, 0]
]

# Initial position of the number to move
x, y = 3, 3

# Move sequence
moves = 'rrldd'

# Function to perform the move
def move_number(matrix, x, y, moves):
    for move in moves:
        if move == 'r':  # Move right
            if y < 4 and matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
                y += 1
        elif move == 'l':  # Move left
            if y > 0 and matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
                y -= 1
        elif move == 'd':  # Move down
            if x < 4 and matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
                x += 1
        elif move == 'u':  # Move up
            if x > 0 and matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
                x -= 1
    return matrix

# Perform the moves
final_matrix = move_number(matrix, x, y, moves)

# Print the final matrix
print(final_matrix)