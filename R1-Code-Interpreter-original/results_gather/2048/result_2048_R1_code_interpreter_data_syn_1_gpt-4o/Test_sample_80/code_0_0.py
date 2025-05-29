# Initial matrix
matrix = [
    [0, 0, 0, 0, 16],
    [4, 32, 0, 32, 4],
    [2, 0, 0, 0, 0],
    [0, 0, 0, 2, 8],
    [0, 0, 0, 8, 0]
]

# Initial position of the number
position = (0, 4)

# Sequence of moves
moves = 'luuddd'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'l':
            if y > 0:
                y -= 1
        elif move == 'r':
            if y < 4:
                y += 1
        elif move == 'u':
            if x > 0:
                x -= 1
        elif move == 'd':
            if x < 4:
                x += 1
        # Check for merging
        if matrix[x][y] == matrix[position[0]][position[1]]:
            matrix[x][y] *= 2
            matrix[position[0]][position[1]] = 0
        else:
            matrix[x][y] = matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
        position = (x, y)
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)