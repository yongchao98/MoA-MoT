# Initial matrix
matrix = [
    [0, 0, 8, 32, 0],
    [8, 0, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [4, 4, 0, 32, 4],
    [2, 32, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 1)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'lrrur'

# Function to perform a move
def move_number(matrix, position, direction):
    row, col = position
    if direction == 'l':
        # Move left
        while col > 0 and matrix[row][col-1] == 0:
            matrix[row][col-1] = matrix[row][col]
            matrix[row][col] = 0
            col -= 1
    elif direction == 'r':
        # Move right
        while col < 4 and matrix[row][col+1] == 0:
            matrix[row][col+1] = matrix[row][col]
            matrix[row][col] = 0
            col += 1
    elif direction == 'u':
        # Move up
        while row > 0 and matrix[row-1][col] == 0:
            matrix[row-1][col] = matrix[row][col]
            matrix[row][col] = 0
            row -= 1
    elif direction == 'd':
        # Move down
        while row < 4 and matrix[row+1][col] == 0:
            matrix[row+1][col] = matrix[row][col]
            matrix[row][col] = 0
            row += 1
    return (row, col)

# Apply the sequence of moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)