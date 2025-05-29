# Initial matrix
matrix = [
    [4, 4],
    [8, 8]
]

# Move sequence
sequence = 'lr'

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

# Function to move right
def move_right(matrix):
    for i in range(2):
        if matrix[i][1] == 0:
            matrix[i][1] = matrix[i][0]
            matrix[i][0] = 0
        elif matrix[i][1] == matrix[i][0]:
            matrix[i][1] *= 2
            matrix[i][0] = 0
    return matrix

# Apply the sequence
for move in sequence:
    if move == 'l':
        matrix = move_left(matrix)
    elif move == 'r':
        matrix = move_right(matrix)

# Print the final matrix
print(matrix)