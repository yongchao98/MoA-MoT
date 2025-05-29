# Initial matrix
matrix = [
    [0, 32, 0, 0],
    [4, 4, 0, 2],
    [16, 0, 0, 32],
    [0, 16, 0, 0]
]

# Move sequence
moves = 'ddll'

# Initial position of the number to move
pos = (2, 3)

# Function to perform the move
def move_number(matrix, pos, direction):
    x, y = pos
    if direction == 'd':
        if x < 3:
            if matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
                x += 1
            elif matrix[x+1][y] == matrix[x][y]:
                matrix[x+1][y] *= 2
                matrix[x][y] = 0
                x += 1
    elif direction == 'l':
        if y > 0:
            if matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
                y -= 1
            elif matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
                y -= 1
    return (x, y)

# Perform the moves
for move in moves:
    pos = move_number(matrix, pos, move)

# Print the final matrix
print(matrix)