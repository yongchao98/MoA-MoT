# Initial matrix
matrix = [
    [8, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 8, 0],
    [0, 0, 0, 0, 0],
    [0, 32, 0, 16, 32]
]

# Initial position of the number to move
x, y = 4, 3

# Sequence of moves
moves = 'uulll'

# Function to perform the move
def move_number(matrix, x, y, moves):
    for move in moves:
        if move == 'u':
            if x > 0:
                x -= 1
        elif move == 'd':
            if x < len(matrix) - 1:
                x += 1
        elif move == 'l':
            if y > 0:
                y -= 1
        elif move == 'r':
            if y < len(matrix[0]) - 1:
                y += 1
    return x, y

# Perform the moves
new_x, new_y = move_number(matrix, x, y, moves)

# Update the matrix
matrix[new_x][new_y] = matrix[4][3]
matrix[4][3] = 0

# Print the updated matrix
print(matrix)