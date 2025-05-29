# Initial matrix
matrix = [
    [0, 0, 8],
    [32, 32, 32],
    [0, 0, 0]
]

# Initial position of the number
x, y = 1, 1

# Move sequence
moves = 'lrl'

# Simulate the moves
for move in moves:
    if move == 'l':
        if y > 0:
            if matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            else:
                matrix[x][y-1], matrix[x][y] = matrix[x][y], matrix[x][y-1]
            y -= 1
    elif move == 'r':
        if y < 2:
            if matrix[x][y+1] == matrix[x][y]:
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            else:
                matrix[x][y+1], matrix[x][y] = matrix[x][y], matrix[x][y+1]
            y += 1

# Print the final matrix
print(matrix)