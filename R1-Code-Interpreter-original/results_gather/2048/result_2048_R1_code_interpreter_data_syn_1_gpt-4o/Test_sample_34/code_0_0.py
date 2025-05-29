# Initial matrix
matrix = [
    [32, 0, 0],
    [2, 4, 0],
    [32, 0, 0]
]

# Move sequence
moves = 'lrl'

# Function to perform a move
def move(matrix, direction):
    if direction == 'l':
        for i in range(3):
            row = [x for x in matrix[i] if x != 0]
            for j in range(len(row) - 1):
                if row[j] == row[j + 1]:
                    row[j] *= 2
                    row[j + 1] = 0
            row = [x for x in row if x != 0]
            matrix[i] = row + [0] * (3 - len(row))
    elif direction == 'r':
        for i in range(3):
            row = [x for x in matrix[i] if x != 0]
            for j in range(len(row) - 1, 0, -1):
                if row[j] == row[j - 1]:
                    row[j] *= 2
                    row[j - 1] = 0
            row = [x for x in row if x != 0]
            matrix[i] = [0] * (3 - len(row)) + row
    return matrix

# Apply the moves
for move_direction in moves:
    matrix = move(matrix, move_direction)

# Print the final matrix
print(matrix)