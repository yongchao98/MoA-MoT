# Initial matrix
matrix = [
    [0, 0, 0, 0, 16],
    [4, 32, 0, 32, 4],
    [2, 0, 0, 0, 0],
    [0, 0, 0, 2, 8],
    [0, 0, 0, 8, 0]
]

# Initial position of the number 16
position = (0, 4)

# Sequence of moves
moves = 'luuddd'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'l':
            if y > 0 and matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
                y -= 1
        elif move == 'u':
            if x > 0 and matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
                x -= 1
        elif move == 'd':
            if x < 4 and matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
                x += 1
        elif move == 'r':
            if y < 4 and matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
                y += 1
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)