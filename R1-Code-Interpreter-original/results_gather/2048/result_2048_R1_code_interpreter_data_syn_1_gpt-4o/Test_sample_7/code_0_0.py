# Initial matrix
matrix = [
    [16, 4],
    [32, 8]
]

# Move sequence
moves = 'rl'

# Function to move right
def move_right(matrix):
    for i in range(2):
        if matrix[i][1] == 0:
            matrix[i][1] = matrix[i][0]
            matrix[i][0] = 0
        elif matrix[i][0] == matrix[i][1]:
            matrix[i][1] *= 2
            matrix[i][0] = 0
    return matrix

# Function to move left
def move_left(matrix):
    for i in range(2):
        if matrix[i][0] == 0:
            matrix[i][0] = matrix[i][1]
            matrix[i][1] = 0
        elif matrix[i][0] == matrix[i][1]:
            matrix[i][0] *= 2
            matrix[i][1] = 0
    return matrix

# Apply moves
for move in moves:
    if move == 'r':
        matrix = move_right(matrix)
    elif move == 'l':
        matrix = move_left(matrix)

# Print the final matrix
print(matrix)